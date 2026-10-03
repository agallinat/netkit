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

<!-- TODO before submitting: add the platforms actually checked. The GitHub
     Actions workflow in .github/workflows/R-CMD-check.yaml covers
     ubuntu-latest (R devel, release, oldrel-1), macOS-latest (release) and
     windows-latest (release), but list a platform here only once a run has
     genuinely passed on it. -->

## Downstream dependencies

There are currently no downstream dependencies for this package.
