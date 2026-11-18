"""Survival by genotype with delayed entry, for any cancer and any data holder.

Phase 1: read a source into the standard tables (delayed_entry.tables), screen them against a study's
settings without reading any outcome (delayed_entry.screen), and lock the filled-in plan (delayed_entry.lock).
Phase 2: the run (delayed_entry.run, with delayed_entry.analysis and the R code in R/), the safety gate
(delayed_entry.gate) and the reports (delayed_entry.reports, delayed_entry.run_report).
Since then: an OMOP CDM adapter, more flat-file formats (SAS, Stata, MAF, VCF), panel sizes from GENIE, the
attribution a source's terms require, and the guided setup (delayed_entry.setup).
"""
__version__ = "0.3.1"
