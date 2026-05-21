# Contributing to SafeInCave

We appreciate contributions of all kinds, whether you are fixing a bug, adding functionality, improving the documentation, or simply opening a discussion. The steps below describe the recommended workflow.

## 1. Start by creating an issue

Please open an issue before you begin working on a change. Use it to explain the problem you want to solve or the feature you would like to add. This helps avoid duplicated work and allows the maintainers to comment on the proposed approach early.

## 2. Fork the repository and clone it locally

Create a fork of the repository on GitHub and then clone your fork to your machine:

```bash
git clone https://github.com/<your-username>/SafeInCave.git
cd SafeInCave
```

Then install the package in editable mode:

```bash
pip install -e .
```

> Note: SafeInCave is intended to run on Linux. If you are using Windows, please work through WSL.

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
