# Gynandromorph storage audit

Allocated regular-file and symlink blocks counted once per inode: 1.399 TiB.
Apparent file sizes: 1.404 TiB. Directory metadata and filesystem overhead are not included.

This is a read-only snapshot. Active run trees and their explicitly traced input dependencies are protected.
Live file sizes are provisional. Historical configuration references are retained as dependency evidence.
An unreferenced intermediate is a review candidate, not permission to remove it.

Use directory_sizes.tsv and largest_files.tsv to identify storage concentration.
cleanup_candidates.tsv contains intermediates, partials and per-run environments, with protection/reference status.
duplicate_candidates.tsv uses size and edge signatures only; it does not prove whole-file identity or reclaimable bytes.
Symlinks and hardlinks are not counted as additional copies of their targets. Removing one hardlink may reclaim no space.

The full compressed inventory and SQLite index remain on Pegasus. Only compact reports are pushed to Git.
No deletion commands or cleanup operations are generated.
