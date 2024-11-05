General Structure
=================

External libraries:

+ the `Eigen <http://eigen.tuxfamily.org/index.php?title=Main_Page>`_ template
  library for linear algebra.  Almost every operation involving matrices and
  vectors is performed through Eigen.  Eigen provides convenient type
  definitions for vectors and matrices (of arbitrary dimensions) and the
  corresponding operations. Have a look
  `here <http://eigen.tuxfamily.org/dox/group__QuickRefPage.html>`_ for a quick
  reference guide to the API and
  at the `getting started guide <http://eigen.tuxfamily.org/dox/GettingStarted.html>`_ to get started.
  Eigen is distributed under the terms of the `Mozilla Public License, v2.0
  <http://opensource.org/licenses/MPL-2.0>`_
+ the `XCFun library <https://xcfun.address.here/missing>`_ by Ulf
  Ekström under the terms of the `GNU General
  Public License, v2.0 <http://opensource.org/licenses/GPL-2.0>`_

 Documentation build:
  
+ the build of this documentation has been copied with minor
  adaptations from the `PCMSolver API
  <https://pcmsolver.link.here>`_. In particular the Sphinx configuration
  script (conf.py) and the Doxygen configuration file (Doxygen.in) and
  the corresponding cmake structure (FindSphinx, find_python_module).

  
