# genie_bpc_crc_synthetic

Made-up data shaped like the GENIE BPC CRC v2.0-public files. No real patients.

Copied unchanged from the CRC BPC delayed-entry project's `fixtures/synthetic_v1`, which that project's
generator made (`crcbpc.synthetic`, 1,500 made-up patients, seed 20261106). The tests here check that the
pipeline reproduces that project's results on these files. Don't edit them by hand.
`truth.csv` is not a BPC file: it holds what only the generator knows, and only the tests read it.
