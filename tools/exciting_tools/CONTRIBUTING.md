# excitingtools Contributing Guide <!-- omit in toc -->

Thank you for investing your time in contributing to our project! Any contribution you make will be reflected in the 
[README](README.md) Contributors list. In this guide you will get an overview of the contribution workflow for opening 
an issue, creating a PR, reviewing, and merging the PR.

## New Contributor Guide

To get an overview of the project, read the [README](README.md). Here are some resources to help you get started with open source contributions:


- [Finding ways to contribute to open source on GitHub](https://docs.github.com/en/get-started/exploring-projects-on-github/finding-ways-to-contribute-to-open-source-on-github)

- [Set up Git](https://docs.github.com/en/get-started/quickstart/set-up-git)

- [GitHub flow](https://docs.github.com/en/get-started/quickstart/github-flow)

- [Collaborating with pull requests](https://docs.github.com/en/github/collaborating-with-pull-requests)


## Contributing to excitingtools

### Create a New Issue

If you spot a problem with the source, tests or docs, check if an issue already exists on Github. In the source code, 
this should be signified by a `TODO(asignee) Number. Description`. If an issue doesn't exist, you can open a new issue.
Please include a clear description of the problem, an indication of which files or features it affects, and a minimum 
working (or failing) example, such that another developer can quickly reproduce the problem. You can either assign 
yourself to the problem (if you wish to work on it), or assign one of the developer team members to it for review.

### Obtaining the Source

1. [Install Git](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git).

2. Fork the repository.

- Using GitHub Desktop:

  - [Getting started with GitHub Desktop](https://docs.github.com/en/desktop/installing-and-configuring-github-desktop/getting-started-with-github-desktop) will guide you through setting up Desktop.

  - Once Desktop is set up, you can use it to [fork the repo](https://docs.github.com/en/desktop/contributing-and-collaborating-using-github-desktop/cloning-and-forking-repositories-from-github-desktop)!

- Using the command line:

  - [Fork the repo](https://docs.github.com/en/github/getting-started-with-github/fork-a-repo#fork-an-example-repository) so that you can make your changes without affecting the original project until you're ready to merge them.

3. Create a working branch and start with your changes:

### Adding a Feature

Each PR should correspond to one feature. This could be a new feature adding functionality to excitingtools, or a fix
for an issue or bug. In any case, please ensure that a branch only contains one feature, such that each merge request 
corresponds to a self-contained change.

Whilst developing, we recommend committing regularly, as all commit messages will be [squashed](https://git-scm.com/book/en/v2/Git-Tools-Rewriting-History)
 into a single commit upon merging. We also recommend pushing changes to your remote regularly. 

Where possible, merge requests should be kept to within 500 lines of code as large submissions become difficult to review.
Furthermore, every new feature or bug fix should be demonstrated with a test. PRs without tests will not be accepted.

### Pull Request

When you're finished with the changes and have written the corresponding test cases, create a [pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request
) targeting the `development` branch of the original repository.

- Be sure to provide clear information in the description indicating:   
  - a). What the feature does,   
  - b). Why it is needed,  
  - c). How to use it.  


- Don't forget to [link the PR to an issue](https://docs.github.com/en/issues/tracking-your-work-with-issues/linking-a-pull-request-to-an-issue) if you are solving one.

- Enable the checkbox to [allow maintainer edits](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/allowing-changes-to-a-pull-request-branch-created-from-a-fork) so the branch can be updated for a merge.

Once you submit your PR, an <span style="font-family:american typewriter; font-size:1em;">**exciting**</span> 
team member will review it:

- It is highly likely that we will ask for changes to be made before a PR can be merged, using pull request comments. 
  We recommend you make these changes in your branch, commit and push.

- As you update your PR and apply changes, mark each conversation as [resolved](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests/commenting-on-a-pull-request#resolving-conversations).

- If you run into any merge issues, checkout this [git tutorial](https://github.com/skills/resolve-merge-conflicts) to help you resolve merge conflicts and other issues.

### Your PR is merged!

Congratulations. Once your PR is merged, your contributions will be added to the project [README](README.md).

## Support

For further support, one can post in the [MatSci exciting forum](https://matsci.org/c/exciting), where one of the
developers will get back to you. One can also contact the lead maintainer, Fabian Peschel via [e-mail](peschelf@physik.hu-berlin.de).
All contributor e-mails are given in [CITATION](CITATION.cff) file.
