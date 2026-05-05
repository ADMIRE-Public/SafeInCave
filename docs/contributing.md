# Contribution Guide

All contributions are welcome: bug fixes, new features, documentation improvements, or just questions. Here's how to get involved.

## 1. Open an issue first

Before starting any work, open an issue to describe what you want to fix or add. This avoids duplicate effort and gives maintainers a chance to give early feedback on the approach.

## 2. Fork the repo and clone your fork

Fork the repository on GitLab, then clone your fork locally:

```bash
git clone https://gitlab.tsn.tno.nl/<your-username>/subsidence-modelling
cd subsidence-modelling
```

Install the project and its dev dependencies:

```bash
uv sync
```

> Note: the project requires Linux. On Windows, use WSL.

## 3. Create a branch

Create a branch in your fork with a short descriptive name:

```bash
git checkout -b my-feature
```

## 4. Make your changes

Run the tests locally to make sure everything passes:

```bash
uv run pytest
```

Run the linter and formatter:

```bash
uv run ruff check .
uv run ruff format .
```

## 5. Open a merge request

Push your branch and open an MR from your fork to the main repository on GitLab. Reference the issue number in the MR description. A maintainer will review and merge it.
