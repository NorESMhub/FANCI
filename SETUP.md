# NorESM / CAM FANCI aerosol — repo setup on Olivia

Setup guide. It covers **two** repos:

| What | Local path | Your fork (`origin`) | Christina's (`christinafork`) |
|---|---|---|---|
| Host model + CIME | `CAM_SEC/` | `<GITHUB_USER>/CAM` | `Trudigard/CAM` |
| FANCI aerosol code (this repo) | `CAM_SEC/src/chemistry/FANCI/` | `<GITHUB_USER>/FANCI` | `Trudigard/FANCI` |

## Prerequisites

- A GitHub account
- Membership of project `nn9560k` on Olivia


---

## 1. Fork on GitHub
Fork both repositories:
- **Host model:** fork  upstream (recommended)
  [`NorESMhub/CAM`](https://github.com/NorESMhub/CAM.git)  OR christinas fork [`Trudigard/CAM`](https://githhttps://github.com/Trudigard/FANCI/blob/saltydust/SETUP.mdub.com/Trudigard/CAM.git)  → `github.com/<GITHUB_USER>/CAM`
- **Aerosol code:** fork [NorESMHub/FANCI](https://github.com/NorESMhub/FANCI)(recommended) OR christinas fork [`Trudigard/FANCI`](https://github.com/Trudigard/FANCI.git) → `github.com/<GITHUB_USER>/FANCI`

## 2. Clone the host model on Olivia

```bash
mkdir -p /cluster/work/projects/nn9560k/$USER
cd /cluster/work/projects/nn9560k/$USER
git clone https://github.com/<GITHUB_USER>/CAM CAM_SEC
cd CAM_SEC
```

## 3. Track Christina's fork and make a working branch

```bash
git remote add christinafork https://github.com/Trudigard/CAM.git
git fetch christinafork
git checkout -b <your_branch_name> christinafork/saltydust
```

> It's recommended to branch with `git checkout -b <name> christinafork/saltydust`. Checking out
> `christinafork/saltydust` directly leaves you in **detached HEAD**, where commits are easy to lose.

## 4. Pull in the external components
[From here]( https://github.com/NorESMhub/noresm3_dev_simulations/wiki/Running-NorESM-on-Olivia)
First
```bash
module purge
module load NRIS/Login
module load Python/3.12.3-GCCcore-13.3.0
module save mod_noresm
```
then next time
```bash
module r mod_noresm

```
Then run:

```bash
./bin/git-fleximod update
```

This checks out `src/chemistry/FANCI` and the other externals at their pinned commits.

## 5. Set up this aerosol repo for development

`git-fleximod` leaves `src/chemistry/FANCI` on a detached commit pointing at Christina's
`FANCI`. Repoint `origin` at your fork so you can commit and push:

```bash
cd src/chemistry/FANCI
git remote set-url origin https://github.com/<GITHUB_USER>/FANCI.git
git remote add christinafork https://github.com/Trudigard/FANCI.git
git fetch christinafork
git checkout -b <your_branch_name> christinafork/saltydust
cd -
```

> After this, `git-fleximod update`
> will warn about the modified external — that is expected; do not let it reset your branch.

## 6. Run the FANCI test suite



From the `CAM_SEC` root:

```bash
./cime/scripts/create_test --xml-category aux_fanci --xml-machine olivia -r /cluster/work/projects/nn9560k/$USER/ -p NN9560K --output-root /cluster/work/projects/nn9560k/$USER/
```

To run a **single** test instead of the whole category, name it explicitly, e.g.:
(COMPSET='SecDust')
```bash
./cime/scripts/create_test SMS_Ln9.ne16pg3_ne16pg3_mtn14.<COMPSET>.olivia_gnu \
  -r /cluster/work/projects/nn9560k/$USER/ -p NN9560K \
  --output-root /cluster/work/projects/nn9560k/$USER/
```

---

## Workflow to contribute code (per now)

**Your own work**

- Commit and push to your fork / your branch — in whichever repo you changed
  (`CAM_SEC/` for host-model or CIME changes, `src/chemistry/FANCI/` for aerosol code).

**Sharing upstream**

- Push your branch to your fork.
- Open a pull request against Christina's repo (`Trudigard/CAM` or `Trudigard/FANCI`).
