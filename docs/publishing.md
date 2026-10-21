# Publishing the prototype repository

The user authorized a **public** GitHub repository named `varquanMS`. The active
GitHub CLI configuration names the account `peijin94`; authentication and account
ownership still need verification before publishing. No GitHub repository was
created by this environment.

Observed local blockers: `api.github.com` cannot be resolved; `gh auth status`
fails; the GitHub connector is disabled by an administrator; and initializing
the existing workspace `.git` fails with `Read-only file system`.
An explicit outside-sandbox `gh repo create` attempt was rejected by the
execution policy because `sandbox_approval=false`, despite the user's permission.

The root project files are ready for Git. The large source MS and tar are under
ignored `data/` paths. A portable Git-history bundle is provided as
`varquanMS-initial.bundle` in this workspace. It contains the initial inspection
commit and prototype implementation, design, tests and measured reports. It
excludes observation payloads, decoded MS/image tables and environment files.

On an authenticated host with GitHub access and a writable destination:

```bash
git clone --branch main /path/to/varquanMS-initial.bundle varquanMS
cd varquanMS
gh auth status
gh repo create varquanMS --public --source=. --remote=origin --push \
  --description 'Adaptive visibility archive experiments for CASA Measurement Sets'
gh repo view --json nameWithOwner,url,visibility
```

Cloning the bundle only restores versioned project files. Copy `data/raw/` and
`data/working/` separately when access to the real observation is needed. Neither
the sample observation nor private environment directories belong in a public
Git push. A software license has not been selected; no license grant is inferred
from repository visibility.
