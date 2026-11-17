"""The script that turns a Markdown document into a Word file. Made-up documents only."""
import importlib.util
import sys
import xml.etree.ElementTree as ET
import zipfile

from helpers import ROOT

spec = importlib.util.spec_from_file_location("md_to_docx", ROOT / "scripts" / "md_to_docx.py")
md = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = md  # its dataclasses look the module up by name
spec.loader.exec_module(md)

W = "{http://schemas.openxmlformats.org/wordprocessingml/2006/main}"
W14 = "{http://schemas.microsoft.com/office/word/2010/wordml}"

SAMPLE = """# A made-up title

Status line with "quotes", it's here, & <odd> text.

## Section

**Proposed:**
- first, with **bold** and *italic*;
  - nested,
    continued here.
- second.

  A paragraph inside the second item.

1. one
2. two
   - under two

   More of two.
3. three

Between the lists.

1. again one

| # | Item | Change |
|---|---------|------|
| 1 | Something | |
| 2 | ☐ Yes ☐ No | |

---

<!-- pagebreak -->

### Last

See [the site](https://example.org/page).

______________________________
"""


def parts(tmp_path, text=SAMPLE):
    src = tmp_path / "doc.md"
    src.write_text(text, encoding="utf-8")
    out = md.convert(src, tmp_path / "doc.docx")
    with zipfile.ZipFile(out) as z:
        return {name: z.read(name) for name in z.namelist()}, out.read_bytes()


def body(tmp_path):
    return ET.fromstring(parts(tmp_path)[0]["word/document.xml"]).find(f"{W}body")


def text_of(element):
    return "".join(t.text or "" for t in element.iter(f"{W}t"))


def test_every_part_is_well_formed_and_listed_and_the_zip_has_no_folder_entries(tmp_path):
    files, _ = parts(tmp_path)
    assert not [n for n in files if n.endswith("/")]
    for name, data in files.items():
        ET.fromstring(data)
    types = files["[Content_Types].xml"].decode()
    for name in files:
        if name.startswith(("word/", "docProps/")) and not name.endswith(".rels"):
            assert f'PartName="/{name}"' in types, name


def test_blocks_in_order_with_list_levels_and_continuations():
    blocks = md.parse(SAMPLE)
    assert [b.kind for b in blocks] == [
        "heading", "para", "heading", "para", "item", "item", "item", "cont",
        "item", "item", "item", "cont", "item", "para", "item", "table", "rule", "pagebreak", "heading",
        "para", "para"]
    nested = blocks[5]
    assert (nested.level, nested.text) == (1, "nested, continued here.")
    assert (blocks[7].level, blocks[7].text) == (0, "A paragraph inside the second item.")
    under_two, more_of_two, three = blocks[10], blocks[11], blocks[12]
    assert (under_two.level, under_two.list_kind) == (1, "bullet")
    assert (more_of_two.kind, more_of_two.level) == ("cont", 0)
    assert (three.list_kind, three.list_id) == ("number", blocks[8].list_id)


def test_each_numbered_list_starts_again_at_one(tmp_path):
    blocks = md.parse(SAMPLE)
    numbered = [b for b in blocks if b.kind == "item" and b.list_kind == "number"]
    assert len({b.list_id for b in numbered}) == 2
    numbering = ET.fromstring(parts(tmp_path)[0]["word/numbering.xml"])
    restarts = numbering.findall(f".//{W}num/{W}lvlOverride/{W}startOverride")
    assert [r.get(f"{W}val") for r in restarts] == ["1", "1"]


def test_table_widths_follow_the_dashes_and_a_form_table_has_rows_to_write_in(tmp_path):
    table = md.parse(SAMPLE)[15]
    widths = md.column_widths(table)
    assert sum(widths) == md.TEXT_W
    assert widths[1] > widths[2] > widths[0]
    tbl = body(tmp_path).find(f"{W}tbl")
    rows = tbl.findall(f"{W}tr")
    assert rows[0].find(f"{W}trPr/{W}tblHeader") is not None
    heights = [r.find(f"{W}trPr/{W}trHeight") for r in rows]
    assert heights[0] is None and all(h is not None for h in heights[1:])


def test_widths_follow_the_text_when_the_dashes_are_equal():
    table = md.parse("| a | b |\n|---|---|\n| x | a much longer cell of text |\n")[0]
    widths = md.column_widths(table)
    assert sum(widths) == md.TEXT_W and widths[1] > widths[0]


def test_each_box_can_be_ticked(tmp_path):
    doc = body(tmp_path)
    assert len(list(doc.iter(f"{W14}checkbox"))) == 2
    assert "Yes" in text_of(doc) and "No" in text_of(doc)


def test_quotes_are_curly_and_odd_characters_survive(tmp_path):
    text = text_of(body(tmp_path))
    assert "“quotes”, it’s here, & <odd> text." in text
    assert '"' not in text and "'" not in text


def test_a_rule_before_a_new_page_is_dropped_and_the_heading_starts_the_page(tmp_path):
    doc = body(tmp_path)
    paragraphs = doc.findall(f"{W}p")
    borders = [p for p in paragraphs if p.find(f"{W}pPr/{W}pBdr") is not None]
    assert borders == []
    last = next(p for p in paragraphs if text_of(p) == "Last")
    assert last.find(f"{W}pPr/{W}pageBreakBefore") is not None


def test_an_item_that_starts_with_a_box_has_no_bullet_and_br_breaks_the_line(tmp_path):
    src = tmp_path / "boxes.md"
    src.write_text("# T\n\n- ☐ First choice\n- ☐ Second choice\n\n| a | b |\n|---|---|\n| ☐ Yes<br>☐ No | |\n",
                   encoding="utf-8")
    md.convert(src, tmp_path / "boxes.docx")
    with zipfile.ZipFile(tmp_path / "boxes.docx") as z:
        doc = ET.fromstring(z.read("word/document.xml")).find(f"{W}body")
    items = [p for p in doc.findall(f"{W}p") if "choice" in text_of(p)]
    assert len(items) == 2 and all(p.find(f"{W}pPr/{W}numPr") is None for p in items)
    cell = doc.find(f"{W}tbl").findall(f"{W}tr")[1].find(f"{W}tc")
    assert len(list(cell.iter(f"{W}br"))) == 1 and len(list(cell.iter(f"{W14}checkbox"))) == 2


def test_a_line_of_underscores_is_text_not_a_rule():
    blocks = md.parse("Above\n\n______________________________\n")
    assert [b.kind for b in blocks] == ["para", "para"]


def test_links_point_outside_and_the_footer_carries_the_title(tmp_path):
    files, _ = parts(tmp_path)
    rels = files["word/_rels/document.xml.rels"].decode()
    assert 'Target="https://example.org/page" TargetMode="External"' in rels
    assert "A made-up title" in files["word/footer1.xml"].decode()


def test_the_same_markdown_gives_the_same_bytes(tmp_path):
    first = parts(tmp_path)[1]
    (tmp_path / "doc.docx").unlink()
    assert parts(tmp_path)[1] == first
