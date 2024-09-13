Kevin T. Chu's Website
============================================================================================

--------------------------------------------------------------------------------------------

This repository contains the source code for my professional website. The website is
implemented in [Jekyll][jekyll].

--------------------------------------------------------------------------------------------

Table of Contents
-----------------

1. [Developer Notes][#1]

   1.1. [Directory Structure][#1.1]

   1.2. [Setting Up for Development][#1.2]

   1.3. [Developing the Website][#1.3]

2. [Deployment Procedures][#2]

3. [References][#3]

--------------------------------------------------------------------------------------------

## 1. Developer Notes

### 1.1. Directory Structure

This project uses the standard Jekyll directory structure with a few additional project
metadata files (e.g., `README.md`).

```
├── README.md          <- this file
├── RELEASE-NOTES.md   <- release notes
├── COPYRIGHT          <- license file
├── VERSION            <- version file
├── Gemfile            <- Ruby package dependency file
├── Gemfile.lock       <- Gem lockfile
├── _config.yml        <- Jekyll configuration
├── _bibliography/     <- directory containing BibTeX file of references
├── _spikes/           <- experimental code
├── _site/             <- Jekyll-generated website (deployable)
├── _tmp/              <- directory for temporary files
├── _*                 <- Jekyll directories containing format, style, data, etc. files
└── all other files    <- website source files
    and directories
```

### 1.2. Setting Up for Development

* Install [Git][git].

* Install [Ruby][ruby].

* Install the Ruby packages required for the project (listed in `Gemfile`).

  ```shell
  $ bundler install
  ```

### 1.3. Developing the Website

#### 1.3.1 Running Local Build and Web Server

* Use the following command to start a Jekyll server that (1) rebuilds the website after
  any source file modifications and (2) runs a local web server for viewing changes.

  ```shell
  $ bundler exec jekyll serve
  ```

#### 1.3.2. Website Data

Data for the website is organized in the `_data` directory as follows.

```
├── _data/
    └── last-updated.yml   <- date that website was last modified
```

#### 1.3.3. Publications

This website uses Jekyll-Scholar to create reference lists and citations for the web pages.
All of the citation data for publications is stored in the `_bibliography/publications.bib`
file. PDF versions of publications are stored in the `research/publications/` directory.

Jekyll-Scholar is configured (in `_config.yml`) to automatically link each citation in
`publications.bib` to a PDF file (when present) named `KEY.pdf` where `KEY` is the key for
the citation. If no PDF file is available for a citation, then no link is generated.

--------------------------------------------------------------------------------------------

## 2. Deployment Procedures

1. Bump the website version in the `VERSION` file.

2. Update the release notes in the `RELEASE-NOTES.md`.

3. Merge changes into the `main` branch on GitHub.

4. The `Jekyll` GitHub Actions workflow will automatically run to build and deploy the
   updated website.

5. (Optional) If the `Jekyll` GitHub Action succeeds, create a new release in GitHub.

--------------------------------------------------------------------------------------------

## 3. References

* [Setting up `Jekyll` GitHub Actions][jekyll-github-actions]

--------------------------------------------------------------------------------------------

[----------------------------------- INTERNAL LINKS -----------------------------------]: #

[#1]: #1-developer-notes
[#1.1]: #11-directory-structure
[#1.2]: #12-setting-up-for-development
[#1.3]: #13-developing-the-website

[#2]: #2-deployment-procedures

[#3]: #3-references

[---------------------------------- REPOSITORY LINKS ----------------------------------]: #

[----------------------------------- EXTERNAL LINKS -----------------------------------]: #

[git]: https://git-scm.com/

[jekyll]: https://jekyllrb.com/

[ruby]: https://www.ruby-lang.org/

[jekyll]: https://jekyllrb.com/

[jekyll-github-actions]: https://jekyllrb.com/docs/continuous-integration/github-actions/
