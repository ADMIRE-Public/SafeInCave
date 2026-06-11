# Contributing to SafeInCave

We appreciate contributions of all kinds, whether you are fixing a bug, adding functionality, improving the documentation, or simply opening a discussion. The steps below describe the recommended workflow.

## 1. Start by creating an issue

Please open an issue before you begin working on a change. Use it to explain the problem you want to solve, or the feature you would like to add, or the documentaion you would like to write. This helps avoid duplicated work and allows the maintainers to comment on the proposed approach early.

## 2. Fork the repository and clone it locally

Create a fork of the repository on GitHub and then clone your fork to your machine.

```bash
git clone https://github.com/<your-username>/SafeInCave.git
cd SafeInCave
```

The user may intend to either (i) develop SafeInCave *source code* or (ii) write *documentation*. Both cases are described below.

### Develop source code

If you want to develop SafeInCave source code, install safeincave conda environment in *development* mode.

```bash
conda env create -f envs/dev.yml
```

Activate envinronment.

```bash
conda activate safeincave-dev
```

Now, all the changes you make in SafeInCave source code should automatically take effect in your simulations.


### Write documentation

If you want to write documentation, install safeincave conda environment in *documentation* mode.

```bash
conda env create -f envs/docs.yml
```

Activate envinronment.

```bash
conda activate safeincave-docs
```

Load documentation.

```bash
mkdocs serve
```

Ctrl+click on the generated URL (usually, http://127.0.0.1:8000/SafeInCave/).


## 3. Create a working branch

Create a new branch in your fork for the changes you plan to make. Choose a short name that clearly describes the purpose of the branch:

```bash
git checkout -b my-feature
```

## 4. Implement and test your changes

Make the required modifications in your branch. Before submitting them, run the test suite locally to check that the existing functionality still works:

```bash
cd tests
python3 -m unittest discover -p "test_*.py"
```

## 5. Submit a pull request

Push your branch to your fork and open a pull request from your fork to the main SafeInCave repository on GitHub. Include the related issue number in the pull request description so the discussion and the proposed changes can be tracked together.

A maintainer will review the contribution and, once everything is ready, merge it into the project.
