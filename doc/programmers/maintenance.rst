Maintenance
===========

Description and how-to for maintenance operations.
Some of the maintenance scripts have been moved to the `pcmsolvermeta
repository <https://gitlab.com/PCMSolver/pcmsolvermeta>`_

Branching Model and Release Process
-----------------------------------

.. warning::
   **Incomplete or outdated information!**

Releases in a ``X.Y.Z`` series are annotated tags on the corresponding branch.

Pull Request Requirements
-------------------------

The project is integrated with `Danger.Systems <http://danger.systems/ruby/>`_.
On each PR, one CI job will run the integration and a `bot <https://github.com/minazobot>`_ will
report which requirements are **not met** in your PR.
These reports can be _warnings_ and _errors_. You will discuss and solve both
of them with the reviewers.
The automatic rules are laid out in the ``Dangerfile`` and are used to enforce an
adequate level of testing, documentation and code quality.

Danger.Systems Warnings
=======================

- PRs classed as Work in Progress.
- Codebase was modified, but no tests were added.
- Nontrivial changes to the codebase, but no documentation added.
- Codebase was modified, but ``CHANGELOG.md`` was not updated.
- Source files were added or removed, but ``.gitattributes`` was not updated.

Danger.Systems Errors
=====================

- Commit message linting, based on some of `these recommendations <https://chris.beams.io/posts/git-commit/>`_:
  - Commit subject is more than one word.
  - Commit subject is no longer than 50 characters.
  - Commit subject and body are separated by an empty line.

- Clean commit history, without merge commits.

- Code style for ``.hpp``, ``.cpp``, ``.h`` files follows the conventions in
  ``.clang-format``.

Bump Version
------------

Version numbering follows the guidelines of `semantic versioning <http://semver.org/>`_
To update, change the relevant field in the ``README.md`` file.

Changelog
---------

We follow the guidelines of `Keep a CHANGELOG <http://keepachangelog.com/>`_
On all **but** the release branches, there is an ``Unreleased`` section
under which new additions should be listed.
To simplify perusal of the ``CHANGELOG.md``, use the following subsections:

1. ``Added`` for new features.
2. ``Changed`` for changes in existing functionality.
3. ``Deprecated`` for once-stable features removed in upcoming releases.
4. ``Removed`` for deprecated features removed in this release.
5. ``Fixed`` for any bug fixes.
6. ``Security`` to invite users to upgrade in case of vulnerabilities.

Updating Eigen Distribution
---------------------------

The C++ linear algebra library Eigen comes bundled with the module. To update
the distributed version one has to:

1. download the desired version of the library to a scratch location. Eigen's
   website is: http://eigen.tuxfamily.org/
2. unpack the downloaded archive;
3. go into the newly created directory and create a build directory;
4. go into the newly created build directory and type the following (remember
   to substitute @PROJECT_SOURCE_DIR@ with the actual path)

   .. code-block:: bash

    cmake .. -DCMAKE_INSTALL_PREFIX=@PROJECT_SOURCE_DIR@/external/eigen3

Remember to commit and push your modifications.

Git Pre-Commit Hooks
--------------------

[Git pre-commit hooks](https://git-scm.com/book/gr/v2/Customizing-Git-Git-Hooks) are used to
keep track of code style and license header in source files.
Code style is checked using ``clang-format``.

.. warning::
   **You need to install ``clang-format`` (v3.9 recommended) to run the code style validation hook!**

License headers are checked using the ``license_maintainer.py`` script and the
header templates for the different languages used in this project.
The Python script checks the ``.gitattributes`` file to determine which license
headers need to be maintained and in which files:

.. code-block:: bash

   src/pedra/pedra_dlapack.F90 !licensefile
   src/solver/*.hpp licensefile=.githooks/LICENSE-C++

The first line specifies that the file in ``src/pedra/pedra_dlapack.F90`` should
not be touched, while the second line states that all ``.hpp`` files in ``src/solver``
should get an header from the template in ``.githooks/LICENSE-C++``
Location of files in ``.gitattributes`` are always specified with respect
to the project root directory.

The hooks are located in the ``.githooks`` subdirectory and **have to be installed by hand**
whenever you clone the repository anew:

.. code-block:: bash

   cd .git/hooks
   cp --symbolic-link ../../.githooks/* .

Installed hooks will **always** be executed. Use ``git commit --no-verify`` to
bypass explicitly the hooks.
