## Submission

This is a new release.

## R CMD check results

0 errors | 0 warnings | 1 note

```
* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Alex Gallinat <alex.gaoc@gmail.com>'
  New submission
```

That note is expected for a first submission.

## Notes on the local check environment

Local `R CMD check --as-cran` additionally reports:

```
* checking for future file timestamps ... NOTE
  unable to verify current time
```

This is not a property of the package. R queries an external time service for
that check, and the service it uses was returning HTTP 502 at the time of
checking. It does not appear when the service is reachable.

## Test environments

* local macOS (aarch64-apple-darwin20), R 4.4.2
* GitHub Actions (`.github/workflows/R-CMD-check.yaml`), all passing:
  * ubuntu-latest, R devel
  * ubuntu-latest, R release
  * ubuntu-latest, R oldrel-1
  * macos-latest, R release
  * windows-latest, R release

The GitHub Actions jobs run `rcmdcheck::rcmdcheck()` with `--as-cran` and
`error_on = "warning"`, so they fail on a warning as well as an error.

## Downstream dependencies

There are currently no downstream dependencies for this package.
