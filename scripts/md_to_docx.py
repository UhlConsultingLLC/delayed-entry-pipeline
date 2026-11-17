"""Turn a Markdown document into a Word file (.docx), with the standard library only.

    python scripts/md_to_docx.py <document.md> [-o <document.docx>]

Written for the study's documents that go to people outside the project, such as the reviewer packet. The
Word file is made from the Markdown, so the Markdown stays the one to edit: change it, then run this again.
The same Markdown always gives the same file, byte for byte.

It reads the part of Markdown those documents use:
- headings (#, ##, ###), paragraphs, bold (**...**), italic (*...*), code (`...`) and links ([text](url));
- bulleted (-, *) and numbered (1.) lists, nested by indentation, and a paragraph that continues a list item
  (indented to the item's text, after a blank line). An item that starts with a box (☐) shows the box
  instead of a bullet;
- <br> for a new line within a paragraph or a table cell;
- pipe tables. When the dashes under the header differ in number from column to column, the columns' widths
  follow them; otherwise they follow the columns' text. A table whose last column is empty in every row is a
  form, so its rows are made tall enough to write in;
- a line of three dashes (a rule), and the comment <!-- pagebreak --> (a new page; Markdown doesn't show it).
  Other one-line comments are left out.

Straight quotes become curly ones, and each ☐ becomes a box that can be ticked in Word. The page is US Letter
with 1-inch margins, and the footer carries the document's first heading and the page number.
"""
from __future__ import annotations

import argparse
import re
import zipfile
from dataclasses import dataclass, field
from pathlib import Path
from xml.sax.saxutils import escape

PAGE_W, PAGE_H, MARGIN = 12240, 15840, 1440  # twentieths of a point: US Letter, 1-inch margins
TEXT_W = PAGE_W - 2 * MARGIN
INDENT = 360  # each list level, a quarter inch
FORM_ROW_HEIGHT = 720  # a form table's rows, at least half an inch
BOX = "☐"
ZIP_DATE = (1980, 1, 1, 0, 0, 0)  # fixed, so the same Markdown gives the same bytes

HEADING = re.compile(r"^(#{1,6}) +(.*?)\s*#*\s*$")
RULE = re.compile(r"^ {0,3}(?:-{3,}|\*{3,})\s*$")
ITEM = re.compile(r"^( *)([-*]|\d+[.)]) +(.*)$")
TABLE_SEP = re.compile(r"^\s*\|?\s*:?-+:?\s*(?:\|\s*:?-+:?\s*)*\|?\s*$")
COMMENT = re.compile(r"^\s*<!--(.*?)-->\s*$")
TOKEN = re.compile(
    r"\*\*(?P<b>.+?)\*\*"
    r"|(?<![\w*])\*(?=\S)(?P<i>[^*]+?)(?<=\S)\*(?![\w*])"
    r"|`(?P<c>[^`]+)`"
    r"|\[(?P<lt>[^\]]+)\]\((?P<lu>[^)\s]+)\)"
    r"|(?P<box>" + BOX + ")"
    r"|(?P<br><br\s*/?>)"
)

W_NS = "http://schemas.openxmlformats.org/wordprocessingml/2006/main"
R_NS = "http://schemas.openxmlformats.org/officeDocument/2006/relationships"
DOC_NS = (f'xmlns:w="{W_NS}" xmlns:r="{R_NS}" '
          'xmlns:w14="http://schemas.microsoft.com/office/word/2010/wordml" '
          'xmlns:mc="http://schemas.openxmlformats.org/markup-compatibility/2006" mc:Ignorable="w14"')
XML_HEAD = '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>\n'
REL_TYPE = "http://schemas.openxmlformats.org/officeDocument/2006/relationships/"


@dataclass
class Block:
    """One block of the document: heading, para, item (a list item), cont (a paragraph continuing a list
    item), table, rule or pagebreak."""
    kind: str
    text: str = ""
    level: int = 0
    list_kind: str = ""  # an item's list: bullet or number
    list_id: int = 0  # numbered lists only: which list, so each starts again at 1
    rows: list = field(default_factory=list)
    dashes: list = field(default_factory=list)


@dataclass
class _Level:
    marker: int  # indent of the list marker
    content: int  # indent of the item's text
    kind: str
    list_id: int


def split_row(line: str) -> list[str]:
    s = line.strip()
    if s.startswith("|"):
        s = s[1:]
    if s.endswith("|") and not s.endswith("\\|"):
        s = s[:-1]
    return [c.strip().replace("\\|", "|") for c in re.split(r"(?<!\\)\|", s)]


def parse(md: str) -> list[Block]:
    """The document's blocks, in order."""
    lines = md.replace("\r\n", "\n").split("\n")
    blocks: list[Block] = []
    stack: list[_Level] = []  # the open list levels
    current: Block | None = None  # the block a following text line continues
    blank = False
    numbered = 0
    i = 0
    while i < len(lines):
        line = lines[i]
        stripped = line.strip()
        indent = len(line) - len(line.lstrip(" "))
        if not stripped:
            current, blank = None, True
            i += 1
            continue
        comment = COMMENT.match(line)
        heading = HEADING.match(line) if indent == 0 else None
        item = ITEM.match(line)
        if comment:
            stack.clear()
            current = None
            if comment.group(1).strip().lower() == "pagebreak":
                blocks.append(Block("pagebreak"))
        elif heading:
            stack.clear()
            current = None
            blocks.append(Block("heading", heading.group(2), level=len(heading.group(1))))
        elif RULE.match(line):
            stack.clear()
            current = None
            blocks.append(Block("rule"))
        elif stripped.startswith("|") and i + 1 < len(lines) and TABLE_SEP.match(lines[i + 1]):
            stack.clear()
            current = None
            header = split_row(line)
            dashes = [len(c.strip(" :")) for c in split_row(lines[i + 1])]
            rows = [header]
            i += 2
            while i < len(lines) and lines[i].strip().startswith("|"):
                row = split_row(lines[i])
                rows.append((row + [""] * len(header))[:len(header)])
                i += 1
            blocks.append(Block("table", rows=rows, dashes=(dashes + [3] * len(header))[:len(header)]))
            blank = False
            continue
        elif item:
            marker, content = len(item.group(1)), item.start(3)
            kind = "number" if item.group(2)[0].isdigit() else "bullet"
            while stack and marker < stack[-1].marker:
                stack.pop()
            if stack and marker < stack[-1].content and stack[-1].kind == kind:
                stack[-1].content = content  # the next item of the same list
            else:
                if stack and marker < stack[-1].content:
                    stack.pop()  # a list of the other kind at the same level ends this one
                stack.append(_Level(marker, content, kind, numbered))
                if kind == "number":
                    numbered += 1
            current = Block("item", item.group(3).strip(), level=len(stack) - 1, list_kind=kind,
                            list_id=stack[-1].list_id)
            blocks.append(current)
        elif current is not None:
            current.text += " " + stripped  # the paragraph or item goes on
        elif stack and blank and indent >= stack[0].content:
            while indent < stack[-1].content:
                stack.pop()
            current = Block("cont", stripped, level=len(stack) - 1)
            blocks.append(current)
        else:
            stack.clear()
            current = Block("para", stripped)
            blocks.append(current)
        blank = False
        i += 1
    return blocks


def smart(text: str, before: str = " ") -> str:
    """Curly quotes for straight ones."""
    out, prev = [], before
    for ch in text:
        if ch == '"':
            ch = "“" if (prev.isspace() or prev in "([{*_/–—") else "”"
        elif ch == "'":
            ch = "‘" if (prev.isspace() or prev in "([{*_“") else "’"
        out.append(ch)
        prev = ch
    return "".join(out)


def plain(text: str) -> str:
    """A line of Markdown as plain text, for the footer and the file's title."""
    text = re.sub(r"\[([^\]]+)\]\([^)]*\)", r"\1", text)
    return smart(text.replace("**", "").replace("`", "").replace("*", ""))


def attr(value: str) -> str:
    return escape(value, {'"': "&quot;"})


class Writer:
    """Makes the runs, and keeps what the document's other parts need: the links and the tick boxes."""

    def __init__(self) -> None:
        self.links: list[str] = []
        self.boxes = 0

    def run(self, text: str, bold=False, italic=False, code=False, style=None) -> str:
        if not text:
            return ""
        props = []
        if style:
            props.append(f'<w:rStyle w:val="{style}"/>')
        if code:
            props.append('<w:rFonts w:ascii="Consolas" w:hAnsi="Consolas" w:cs="Consolas"/>')
        if bold:
            props.append("<w:b/><w:bCs/>")
        if italic:
            props.append("<w:i/><w:iCs/>")
        rpr = f"<w:rPr>{''.join(props)}</w:rPr>" if props else ""
        return f'<w:r>{rpr}<w:t xml:space="preserve">{escape(text)}</w:t></w:r>'

    def link(self, text: str, url: str, bold: bool, italic: bool) -> str:
        self.links.append(url)
        runs = self.run(smart(text), bold, italic, style="Hyperlink")
        return f'<w:hyperlink r:id="rIdLink{len(self.links)}" w:history="1">{runs}</w:hyperlink>'

    def checkbox(self) -> str:
        self.boxes += 1
        font = '<w:rFonts w:ascii="MS Gothic" w:eastAsia="MS Gothic" w:hAnsi="MS Gothic" w:hint="eastAsia"/>'
        return ('<w:sdt><w:sdtPr>'
                f'<w:rPr>{font}</w:rPr><w:id w:val="{self.boxes}"/>'
                '<w14:checkbox><w14:checked w14:val="0"/>'
                '<w14:checkedState w14:val="2612" w14:font="MS Gothic"/>'
                '<w14:uncheckedState w14:val="2610" w14:font="MS Gothic"/></w14:checkbox>'
                f'</w:sdtPr><w:sdtContent><w:r><w:rPr>{font}</w:rPr><w:t>{BOX}</w:t></w:r></w:sdtContent></w:sdt>')

    def inline(self, text: str, bold=False, italic=False, before=" ") -> str:
        """The runs for a line of Markdown text."""
        out, pos = [], 0
        for m in TOKEN.finditer(text):
            if m.start() > pos:
                out.append(self.run(smart(text[pos:m.start()], text[pos - 1] if pos else before), bold, italic))
            if m.group("b") is not None:
                out.append(self.inline(m.group("b"), True, italic, "*"))
            elif m.group("i") is not None:
                out.append(self.inline(m.group("i"), bold, True, "*"))
            elif m.group("c") is not None:
                out.append(self.run(m.group("c"), bold, italic, code=True))
            elif m.group("lt") is not None:
                out.append(self.link(m.group("lt"), m.group("lu"), bold, italic))
            elif m.group("br") is not None:
                out.append("<w:r><w:br/></w:r>")
            else:
                out.append(self.checkbox())
            pos = m.end()
        if pos < len(text):
            out.append(self.run(smart(text[pos:], text[pos - 1] if pos else before), bold, italic))
        return "".join(out)


def paragraph(inner: str, style: str | None = None, *, num: int | None = None, ilvl: int = 0,
              left: int | None = None, keep_next=False, page_break=False, before: int | None = None,
              after: int | None = None, border=False) -> str:
    ppr = []
    if style:
        ppr.append(f'<w:pStyle w:val="{style}"/>')
    if keep_next:
        ppr.append("<w:keepNext/>")
    if page_break:
        ppr.append("<w:pageBreakBefore/>")
    if num is not None:
        ppr.append(f'<w:numPr><w:ilvl w:val="{ilvl}"/><w:numId w:val="{num}"/></w:numPr>')
    if border:
        ppr.append('<w:pBdr><w:bottom w:val="single" w:sz="6" w:space="1" w:color="BFBFBF"/></w:pBdr>')
    spacing = " ".join(f'w:{k}="{v}"' for k, v in (("before", before), ("after", after)) if v is not None)
    if spacing:
        ppr.append(f"<w:spacing {spacing}/>")
    if left is not None:
        ppr.append(f'<w:ind w:left="{left}"/>')
    props = f"<w:pPr>{''.join(ppr)}</w:pPr>" if ppr else ""
    return f"<w:p>{props}{inner}</w:p>"


def column_widths(block: Block) -> list[int]:
    """Each column's width: by the dashes under the header when they differ, otherwise by the text."""
    n = len(block.rows[0])
    if len(set(block.dashes)) > 1:
        weights = list(block.dashes)
    else:
        weights = [min(max(max(len(r[j]) for r in block.rows), 4), 50) for j in range(n)]
    total = sum(weights)
    widths = [TEXT_W * w // total for w in weights]
    widths[-1] += TEXT_W - sum(widths)
    return widths


def table(block: Block, w: Writer) -> str:
    widths = column_widths(block)
    body = block.rows[1:]
    form = bool(body) and all(not row[-1] for row in body)
    borders = "".join(f'<w:{side} w:val="single" w:sz="4" w:space="0" w:color="BFBFBF"/>'
                      for side in ("top", "left", "bottom", "right", "insideH", "insideV"))
    out = [f'<w:tbl><w:tblPr><w:tblW w:w="{TEXT_W}" w:type="dxa"/><w:tblBorders>{borders}</w:tblBorders>'
           '<w:tblLayout w:type="fixed"/><w:tblCellMar><w:left w:w="100" w:type="dxa"/>'
           '<w:right w:w="100" w:type="dxa"/></w:tblCellMar></w:tblPr>',
           "<w:tblGrid>" + "".join(f'<w:gridCol w:w="{x}"/>' for x in widths) + "</w:tblGrid>"]
    for r, row in enumerate(block.rows):
        head = r == 0
        trpr = "<w:cantSplit/>" + ("<w:tblHeader/>" if head else "")
        if form and not head:
            trpr += f'<w:trHeight w:val="{FORM_ROW_HEIGHT}" w:hRule="atLeast"/>'
        out.append(f"<w:tr><w:trPr>{trpr}</w:trPr>")
        for j, text in enumerate(row):
            shade = '<w:shd w:val="clear" w:color="auto" w:fill="F2F2F2"/>' if head else ""
            out.append(f'<w:tc><w:tcPr><w:tcW w:w="{widths[j]}" w:type="dxa"/>{shade}</w:tcPr>'
                       f"{paragraph(w.inline(text, bold=head), 'TableText')}</w:tc>")
        out.append("</w:tr>")
    out.append("</w:tbl>")
    return "".join(out)


def document_body(blocks: list[Block], w: Writer) -> str:
    out = []
    page_break = False
    for k, b in enumerate(blocks):
        prev = blocks[k - 1] if k else None
        nxt = blocks[k + 1] if k + 1 < len(blocks) else None
        if b.kind == "pagebreak":
            page_break = True
            continue
        if b.kind == "rule" and nxt is not None and nxt.kind == "pagebreak":
            continue  # the new page is division enough
        after_table = prev is not None and prev.kind == "table"
        if b.kind == "heading":
            out.append(paragraph(w.inline(b.text), f"Heading{min(b.level, 3)}", page_break=page_break))
        elif b.kind == "para":
            out.append(paragraph(w.inline(b.text), keep_next=nxt is not None and nxt.kind in ("item", "table"),
                                 page_break=page_break, before=160 if after_table else None))
        elif b.kind in ("item", "cont"):
            last = nxt is None or nxt.kind not in ("item", "cont")
            after = 160 if last else 40
            if b.kind == "item" and b.text.startswith(BOX):
                out.append(paragraph(w.inline(b.text), "ListParagraph", left=INDENT * b.level,
                                     page_break=page_break, after=after))  # the box stands for the bullet
            elif b.kind == "item":
                num = 1 if b.list_kind == "bullet" else 2 + b.list_id
                out.append(paragraph(w.inline(b.text), "ListParagraph", num=num, ilvl=min(b.level, 2),
                                     page_break=page_break, after=after))
            else:
                out.append(paragraph(w.inline(b.text), "ListParagraph", left=INDENT * (b.level + 1),
                                     page_break=page_break, before=80, after=after))
        elif b.kind == "table":
            if page_break:
                out.append(paragraph("", page_break=True))
            out.append(table(b, w))
        elif b.kind == "rule":
            out.append(paragraph("", border=True, page_break=page_break, before=0, after=120))
        page_break = False
    if blocks and blocks[-1].kind == "table":
        out.append(paragraph(""))
    return "".join(out)


def styles() -> str:
    def heading(n: int, size: int, before: int, after: int) -> str:
        return (f'<w:style w:type="paragraph" w:styleId="Heading{n}"><w:name w:val="heading {n}"/>'
                '<w:basedOn w:val="Normal"/><w:next w:val="Normal"/><w:uiPriority w:val="9"/><w:qFormat/>'
                f'<w:pPr><w:keepNext/><w:keepLines/><w:spacing w:before="{before}" w:after="{after}"/>'
                f'<w:outlineLvl w:val="{n - 1}"/></w:pPr>'
                f'<w:rPr><w:b/><w:bCs/><w:color w:val="1F3864"/><w:sz w:val="{size}"/><w:szCs w:val="{size}"/>'
                '</w:rPr></w:style>')
    return (XML_HEAD + f'<w:styles xmlns:w="{W_NS}">'
            '<w:docDefaults><w:rPrDefault><w:rPr>'
            '<w:rFonts w:ascii="Calibri" w:eastAsia="Calibri" w:hAnsi="Calibri" w:cs="Calibri"/>'
            '<w:sz w:val="22"/><w:szCs w:val="22"/><w:lang w:val="en-US"/></w:rPr></w:rPrDefault>'
            '<w:pPrDefault><w:pPr><w:spacing w:after="120" w:line="264" w:lineRule="auto"/></w:pPr>'
            '</w:pPrDefault></w:docDefaults>'
            '<w:style w:type="paragraph" w:default="1" w:styleId="Normal"><w:name w:val="Normal"/><w:qFormat/>'
            '</w:style>'
            + heading(1, 36, 0, 160) + heading(2, 28, 240, 120) + heading(3, 24, 240, 80) +
            '<w:style w:type="paragraph" w:styleId="ListParagraph"><w:name w:val="List Paragraph"/>'
            '<w:basedOn w:val="Normal"/><w:uiPriority w:val="34"/><w:qFormat/>'
            '<w:pPr><w:spacing w:after="40"/></w:pPr></w:style>'
            '<w:style w:type="paragraph" w:customStyle="1" w:styleId="TableText"><w:name w:val="Table Text"/>'
            '<w:basedOn w:val="Normal"/><w:qFormat/><w:pPr><w:spacing w:before="40" w:after="40" w:line="240" '
            'w:lineRule="auto"/></w:pPr><w:rPr><w:sz w:val="20"/><w:szCs w:val="20"/></w:rPr></w:style>'
            '<w:style w:type="paragraph" w:styleId="Footer"><w:name w:val="footer"/><w:basedOn w:val="Normal"/>'
            f'<w:pPr><w:tabs><w:tab w:val="right" w:pos="{TEXT_W}"/></w:tabs><w:spacing w:after="0"/></w:pPr>'
            '<w:rPr><w:color w:val="595959"/><w:sz w:val="18"/><w:szCs w:val="18"/></w:rPr></w:style>'
            '<w:style w:type="character" w:styleId="Hyperlink"><w:name w:val="Hyperlink"/>'
            '<w:rPr><w:color w:val="0563C1"/><w:u w:val="single"/></w:rPr></w:style>'
            '</w:styles>')


def numbering(n_numbered: int) -> str:
    def level(i: int, fmt: str, text: str) -> str:
        return (f'<w:lvl w:ilvl="{i}"><w:start w:val="1"/><w:numFmt w:val="{fmt}"/><w:lvlText w:val="{text}"/>'
                f'<w:lvlJc w:val="left"/><w:pPr><w:ind w:left="{INDENT * (i + 1)}" w:hanging="{INDENT}"/>'
                '</w:pPr></w:lvl>')
    bullets = "".join(level(i, "bullet", ch) for i, ch in enumerate("•–◦"))
    numbers = "".join(level(i, "decimal", f"%{i + 1}.") for i in range(3))
    nums = '<w:num w:numId="1"><w:abstractNumId w:val="0"/></w:num>' + "".join(
        f'<w:num w:numId="{2 + k}"><w:abstractNumId w:val="1"/>'
        '<w:lvlOverride w:ilvl="0"><w:startOverride w:val="1"/></w:lvlOverride></w:num>'
        for k in range(n_numbered))
    return (XML_HEAD + f'<w:numbering xmlns:w="{W_NS}">'
            f'<w:abstractNum w:abstractNumId="0"><w:multiLevelType w:val="hybridMultilevel"/>{bullets}'
            '</w:abstractNum>'
            f'<w:abstractNum w:abstractNumId="1"><w:multiLevelType w:val="hybridMultilevel"/>{numbers}'
            f'</w:abstractNum>{nums}</w:numbering>')


def footer(title: str) -> str:
    def field_(code: str) -> str:
        return f'<w:fldSimple w:instr=" {code} "><w:r><w:t>1</w:t></w:r></w:fldSimple>'
    return (XML_HEAD + f'<w:ftr xmlns:w="{W_NS}" xmlns:r="{R_NS}"><w:p><w:pPr><w:pStyle w:val="Footer"/></w:pPr>'
            f'<w:r><w:t xml:space="preserve">{escape(title)}</w:t></w:r><w:r><w:tab/></w:r>'
            f'<w:r><w:t xml:space="preserve">Page </w:t></w:r>{field_("PAGE")}'
            f'<w:r><w:t xml:space="preserve"> of </w:t></w:r>{field_("NUMPAGES")}</w:p></w:ftr>')


def build(md: str) -> dict[str, str]:
    """The Word file's parts, by name."""
    blocks = parse(md)
    title = next((plain(b.text) for b in blocks if b.kind == "heading" and b.level == 1), "")
    w = Writer()
    body = document_body(blocks, w)
    n_numbered = len({b.list_id for b in blocks if b.kind == "item" and b.list_kind == "number"})
    section = ('<w:sectPr><w:footerReference w:type="default" r:id="rIdFooter"/>'
               f'<w:pgSz w:w="{PAGE_W}" w:h="{PAGE_H}"/>'
               f'<w:pgMar w:top="{MARGIN}" w:right="{MARGIN}" w:bottom="{MARGIN}" w:left="{MARGIN}" '
               'w:header="720" w:footer="720" w:gutter="0"/></w:sectPr>')
    links = "".join(f'<Relationship Id="rIdLink{k}" Type="{REL_TYPE}hyperlink" Target="{attr(url)}" '
                    'TargetMode="External"/>' for k, url in enumerate(w.links, 1))
    word = "application/vnd.openxmlformats-officedocument.wordprocessingml"
    return {
        "[Content_Types].xml": XML_HEAD +
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        f'<Override PartName="/word/document.xml" ContentType="{word}.document.main+xml"/>'
        f'<Override PartName="/word/styles.xml" ContentType="{word}.styles+xml"/>'
        f'<Override PartName="/word/numbering.xml" ContentType="{word}.numbering+xml"/>'
        f'<Override PartName="/word/settings.xml" ContentType="{word}.settings+xml"/>'
        f'<Override PartName="/word/footer1.xml" ContentType="{word}.footer+xml"/>'
        '<Override PartName="/docProps/core.xml" '
        'ContentType="application/vnd.openxmlformats-package.core-properties+xml"/></Types>',
        "_rels/.rels": XML_HEAD +
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        f'<Relationship Id="rId1" Type="{REL_TYPE}officeDocument" Target="word/document.xml"/>'
        '<Relationship Id="rId2" '
        'Type="http://schemas.openxmlformats.org/package/2006/relationships/metadata/core-properties" '
        'Target="docProps/core.xml"/></Relationships>',
        "docProps/core.xml": XML_HEAD +
        '<cp:coreProperties xmlns:cp="http://schemas.openxmlformats.org/package/2006/metadata/core-properties" '
        'xmlns:dc="http://purl.org/dc/elements/1.1/" xmlns:dcterms="http://purl.org/dc/terms/" '
        'xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">'
        f'<dc:title>{escape(title)}</dc:title></cp:coreProperties>',
        "word/_rels/document.xml.rels": XML_HEAD +
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        f'<Relationship Id="rIdStyles" Type="{REL_TYPE}styles" Target="styles.xml"/>'
        f'<Relationship Id="rIdNumbering" Type="{REL_TYPE}numbering" Target="numbering.xml"/>'
        f'<Relationship Id="rIdSettings" Type="{REL_TYPE}settings" Target="settings.xml"/>'
        f'<Relationship Id="rIdFooter" Type="{REL_TYPE}footer" Target="footer1.xml"/>'
        f'{links}</Relationships>',
        "word/document.xml": XML_HEAD + f'<w:document {DOC_NS}><w:body>{body}{section}</w:body></w:document>',
        "word/styles.xml": styles(),
        "word/numbering.xml": numbering(n_numbered),
        "word/settings.xml": XML_HEAD +
        f'<w:settings xmlns:w="{W_NS}"><w:defaultTabStop w:val="720"/>'
        '<w:characterSpacingControl w:val="doNotCompress"/><w:compat>'
        '<w:compatSetting w:name="compatibilityMode" w:uri="http://schemas.microsoft.com/office/word" w:val="15"/>'
        '</w:compat></w:settings>',
        "word/footer1.xml": footer(title),
    }


def convert(markdown: Path, out: Path) -> Path:
    parts = build(markdown.read_text(encoding="utf-8"))
    with zipfile.ZipFile(out, "w", zipfile.ZIP_DEFLATED) as z:
        for name, text in parts.items():
            info = zipfile.ZipInfo(name, date_time=ZIP_DATE)
            info.compress_type = zipfile.ZIP_DEFLATED
            z.writestr(info, text.encode("utf-8"))
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description="Turn a Markdown document into a Word file (.docx).")
    ap.add_argument("markdown", type=Path)
    ap.add_argument("-o", "--output", type=Path, help="the Word file (default: the Markdown's name, .docx)")
    args = ap.parse_args()
    out = convert(args.markdown, args.output or args.markdown.with_suffix(".docx"))
    print(f"Wrote {out}")


if __name__ == "__main__":
    main()
