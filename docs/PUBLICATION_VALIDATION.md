# Publication validation and recovery review

Review date: 14 October 2026. Recovery branch: pegasus-publication-recovery-20261014T183358Z.

Syntax and selected dependency imports passed for the original validation checkout. No unit suite was available there. These checks do not validate the recovered historical analysis scripts or reproduce figures.

The recovery contains 301 candidate files: three from the Git checkout and 298 from the analysis workspace. The checkout candidates are a hard-coded joint-site-count audit and two Figure 2 submission wrappers. The latter reference the six B/C/D/E/G/H panels; compare their scheduler settings with current maintained wrappers before replacing anything. The remaining candidates include historical optix cohort/filter revisions, figure alternatives, NDJ modules and tests, and scripts inside analysis directories. A filename or recent modification is not evidence that a version produced an accepted figure.

Keep these candidates on the recovery branch until each retained workflow has an accepted input/reference manifest, dependency closure, and figure/output comparison. Do not bulk merge historical filters into the current joint-recall workflow. No candidate has been promoted by this review.

Publication scope: maintained source, configuration, dependency definitions, reproduction commands, and small intentional test fixtures. Raw scientific data, intermediates, caches, installed environments, and generated run directories stay outside the code release. Passing syntax/tests is distinct from reproducing accepted scientific outputs.
