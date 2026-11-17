# Security policy

## Reporting a problem

Please report a security problem privately: on this repository's **Security** tab, choose **Report a
vulnerability**. Don't open a public issue for it.

## Never post data

This pipeline is built for patient-level clinical and genomic data held under data-use agreements. In an issue,
a pull request or any other public place, never post:
- patient-level data, or any value from a data row;
- a count smaller than the data-use agreement allows;
- credentials or access tokens, such as a Synapse personal access token;
- a local path or file name from a machine that holds real data.

Describe a problem with made-up data, as the tests do.

## What this repository holds

Code, settings, documentation and made-up test data: no real patient data. `tests/test_repository_guard.py`
runs on every push and fails if a committed file:
- looks like a data file outside the made-up fixtures (`tests/fixtures/`);
- is over 1 MB;
- holds a local user folder's path;
- looks like a credential.

## Supported versions

Only the latest version on `main` is maintained.
