# Publishing the prototype repository

The public repository is [peijin94/vamp](https://github.com/peijin94/vamp).
It was created and pushed on 2026-10-22 after outside-sandbox GitHub access
became available. Public visibility and the initial successful Actions run were
verified. Earlier network/approval failures are resolved for publication.

The project directory now has a normal Git checkout tracking `origin/main`.
It was initialized from the published history without overwriting working files.
Git mutations may still require outside-sandbox execution because `.git` is
protected inside the restricted sandbox. The original portable bundle remains a
historical snapshot of the initial two commits, not the latest branch.

The root project files are ready for Git. The large source MS and tar are under
ignored `data/` paths. A portable Git-history bundle is provided as
`vamp-initial.bundle` in this workspace. It contains the initial inspection
commit and prototype implementation, design, tests and measured reports. It
excludes observation payloads, decoded MS/image tables and environment files.

Clone the current repository:

```bash
git clone https://github.com/peijin94/vamp.git
```

Cloning the bundle only restores versioned project files. Copy `data/raw/` and
`data/working/` separately when access to the real observation is needed. Neither
the sample observation nor private environment directories belong in a public
Git push. A software license has not been selected; no license grant is inferred
from repository visibility.
