.. _sec-contributing-4:

****************************
Guidelines for Contributors
****************************

We encourage contributions to the ShakeMap project. Before
you begin, however, we recommend you contact us (either
directly or, preferrably, through the *issues* feature on
GitHub) to discuss your plans so that we can try to
coordinate activities among contributors and ourselves.

See our LICENSE.md file for details of ShakeMap licensing.
Any contributions will fall under the terms of that license,
and by submitting a pull request you are agreeing to
place your code under our license.
In particular, it is a public domain license, which means
you effectively give up all rights to your contributions,
including to someone who uses or sells the code commercially.

Unit tests are required for all contributed code, and must
provide significant coverage (>90%). All existing tests must
pass or the PR will be rejected.

Documentation is also required. In addition to being 
well-commented, make sure your code has
`Google-style docstrings
<http://sphinxcontrib-napoleon.readthedocs.io/en/latest/example_google.html>`_
describe functions and all inputs, 
outputs, and possible exceptions raised. See the section
"Coding Style" below for more.

There is a program, *makedocs.py* at the top level of the
repository. It will rebuild the documentation. It will also
change or create many files. Make sure to run the program 
(and review the documentation) before submitting your pull
request.

See below for details on:

- Coding style
- Development workflow
- Coremods

Coding Style
==========================

In general, follow the `PEP-8 style guide
<https://www.python.org/dev/peps/pep-0008/>`_.
Modifications/details/re-emphasis:

**General**

- Put spaces after commas
- *Seriously*, put spaces after commas
- Put spaces after the beginning '#' of comments
- Put spaces around mathematical operators unless they are part of the
  argument list of a function.
- Comment heavily with useful comments
- Document everything
- Keep imports on separate lines to the extent that this is reasonable
- `autopep8 <https://github.com/hhatto/autopep8>`_ is useful for fixing
  simple whitespace issues (we do not use the more aggressive options).
- Some of us also use the `spyder
  <https://pythonhosted.org/spyder/>`_ IDE, which helps catch many errors
  and style issues. It also has an option for automatically removing
  trailing whitespace. But there are also Python style/syntax checking
  modules for many editors, including vim and emacs.
- Two blank lines between functions. One blank line between methods of a 
  class.
- Put spaces after commas

**Names**

Program names should be:

- Lowercase, underscored when necessary
- Without .py suffix (``model``, not ``model.py``)
- Verbs (``model``, not ``modeler``)

Package names should be:

- lowercased, underscored when necessary
- nouns
- reasonably descriptive

Module names should be:

- Lowercased, underscored when necessary
- Nouns
- Reasonably descriptive

Class names should be:

- Camel-cased, starting with uppercase character (``MapCalculator``)
- Nouns

Instance method names should be:

- Camel-cased, starting with lowercase character [calcRxDistance()]
  (Unless superclass uses another style. In that case, follow the 
  superclass style.)
- Verbs
- Reasonably descriptive (``calcDistance`` is ok, does not need to 
  be ``calculateDistance``)
- If not meant to be used by external users, hidden with single 
  or double underscores

Class method names should be:

- Camel-cased, starting with lowercase character (``calcRxDistance``)
- If not meant to be used by external users, hidden with single or 
  double underscores.

Module or program function names should be:

- Lowercased, underscored when necessary
- If not meant to be used by external users, hidden with single or 
  double underscores.

Variable names should be:

- Lowercased, underscored when necessary
- Descriptive when not used as iterators
- Nouns

Global variable/constant names should be:

- Used sparingly, if at all (especially variables)
- ALL_CAPS with underscores if necessary
- Hidden with leading single or double underscores if not meant to 
  be used outside their module (which, really, they shouldn’t be.) 
  The ALL_CAPS convention should be enough to tell people "do not 
  modify this."

Class variable names should be:

- ALL_CAPS with underscores if necessary
- If not meant to be public, probably hidden with single or double 
  underscores and wrapped with getter function (setter functions 
  are probably bad here).  The ALL_CAPS convention should be enough 
  to tell people "do not modify this."

Class instance variable names should be:

- Lowercased, underscored when necessary (fault_distance)
- If meant to be public, probably hidden with single or double 
  underscores and wrapped with getter function (setter functions 
  are probably bad here).
- Possibly hidden with single or double underscores anyway just 
  to keep them private.

**Misc**

- Please follow the 80 character line limit when possible
- 4 *space* indentation for code blocks (not tabs)--most Python 
  aware editors should support this out of the box
- Whitespace and indentation should follow PEP 8 style guide
- Single letter class and function names are discouraged; 
  single-letter variable names should be restricted to things 
  like indices in loops and arrays and should never be "l", 
  "I", or "O" (lowercase L, uppercase i, and uppercase o, 
  respectively). 

**Exception Handling**

- Prefer the built-in Python exceptions where applicable.  The 
  full list is found `here
  <https://docs.python.org/3.5/library/exceptions.html#exceptions.IOError>`_, 
  but the exceptions most likely to be used are:

  - KeyError -- Invalid key (as in dict or dict-like object)
  - IndexError -- Invalid index (asking for 11th element of 10 element sequence)
  - TypeError -- Operation or function applied to wrong type.
  - ValueError -- Value out of range (i.e., magnitude > 10)
  - FileNotFoundError -- Kinda self-explanatory

- Regardless of Exception type, use a descriptive error message 
  so the user or developer has a good idea of what exactly the 
  problem was.

**File layout**

- Imports:

  - Stdlib imports first
  - Third-party imports second
  - Local imports last

- All globals should be declared at the top of the file 
  after the imports

- Following the globals, the primary class and/or functions 
  of a module should come first. Secondary and helper classes 
  should follow the primary class. Helper functions and other 
  "invisible" stuff should follow. The ``"__main__"`` block (if 
  any) should come last

**Example**

Below is a simple example module with a module-level method
(``convenience_calculator``) and a class (``DummyClass``)
with a class method (``doThingsWithMagnitude``). For real code
examples, see
`the coremods <https://github.com/usgs/shakemap/tree/master/shakemap/coremods)>`_
directory on Github.

.. code-block:: python

  #!/usr/bin/env python

  CONSTANT = 5.0


  def convenience_calculator(value):
      """
      Return the input value multiplied by 5.0.
    
      Args:
          value (float): Input numeric value.

      Returns:
          float: The product of the input value and the number 5.0.
      """
      return value * CONSTANT


  class DummyClass(object):
      """
      This is the one-line description of this class.

      This is the multi-line,
      detailed description of the purpose of this class.
      """
      CALC_VERSION = '1.1'
      __HIDDEN_CLASS_VARIABLE = 2
      def __init__(self, calc_string='not zero'):
          """
          Create an instance of the DummyClass class.

          Args:
              calc_string (str): A string.
          """
          self.calc_string = calc_string
  
      def doThingsWithMagnitude(self, mag, mystr):
          """
          Do mysterious things with magnitude.
  
          More detailed description of the weird and 
          wonderful things that will be done with magnitude.
             
          Args:
              mag (float): Input numeric magnitude value, should 
                  be 0 < mag < 10.
              mystr (str): String input that isn't used.

          Returns:
              float: 0.0 or (CONSTANT * mag), depending on the calc_string
              defined in the constructor.

          Raises: 
              ValueError: If input magnitude is outside the 
                  accepted range.
          """
          if mag <= 0 or mag >= 10:
              msg = 'Input magnitudes must be between 0 and 10.'
              raise ValueError(msg)
          if self.calc_string == 'not zero':
              calc_result = mag * CONSTANT
          else:
              calc_result = 0.0
          return calc_result

Workflow
=========

Below is a description of our Git workflow. This workflow
is an adaption of a fairly common set of procedures for working with
GitHub. :num:`Figure #workflow-figure` illustrates the process described
below.

.. _workflow-figure:

.. figure:: _static/workflow.*
   :width: 710
   :align: center

   ShakeMap Git workflow.

If they have not previously worked on this project, the developer 
must first fork the main USGS ShakeMap source code repository. This 
fork is 
later called the remote **origin**. The developer then clones the fork 
to the developer’s local development workstation. By cloning from a 
remote repository on GitHub, the local repository automatically sets 
up an origin remote reference. It is a good idea to manually define 
an **upstream** remote reference at this time as well::

  $ git clone git@github.com:username/shakemap.git
  $ cd shakemap
  $ git remote add upstream https://github.com/usgs/shakemap.git

Working in their local repository, the developer creates a feature 
branch based off the master branch and begins work. Source code is 
modified and incremental commits are made against the local feature 
branch::

  $ git checkout master
  $ git branch branch-name
  $ git checkout branch-name
  $ vim file1 file2 file3
  $ git commit -am 'Modified three files to implement ticket-number.'
  $ vim test1 test2 test3
  $ git commit -am 'Wrote/updated tests for files changed.'

It is important to note that completing the feature involves both 
completing and testing the feature implementation. Having an automated 
test framework helps reduce regression tests moving forward. We use
Travis CI and tests are automatically initiated when a pull request is
made. Unit tests for new code are required or the developer's pull
request will be rejected. Tests and test data are in the *tests* 
sub-directory. Running tests before submitting a pull request is
strongly recommended::

  $ py.test --cov=. --cov-report html

Assuming the tests pass, this command will place a coverage report in
*htmlcov/index.html* that links to the covered files and shows which 
lines are covered and which are not. Please strive for complete 
coverage -- our goal is to keep project coverage above 90 percent.

When the developer completes work for this feature, they first integrate 
any changes contributed by other developers (i.e., changes that were made
to the **upstream** master while they were working on their branch), and then 
push the feature branch 
back to the remote origin. This is done by pulling changes in the **upstream** 
master branch down to the local master, and then rebasing the local feature 
branch against the new local master branch::

  $ git checkout master
  $ git pull --rebase upstream master
  $ git checkout branch-name
  $ git rebase master
  $ git push origin branch-name

Obviously, if there are conflicts while rebasing, the developer must
resolve them before proceeding. Conflicts can be minimized by 
communication with the ShakeMap team through the issues feature of 
GitHub.

The developer now creates a pull request for this feature. This is done by 
logging into their account on GitHub, navigating to their fork of the 
repository, and clicking on the "Pull Request" button on that page. Once
the pull request has been merged by the ShakeMap repositiory admins,
the developer may delete the feature branch on GitHub and their local 
maching. Then, the developer should rebase their master branch from 
the **upstream** master (which now contains their merged feature)
and push them to their **origin** master::

  $ git pull --rebase upstream master
  $ git push origin master

The developer's repository is then ready to begin work on a new
branch.

For the uninitiated this workflow may seem a bit convoluted, but it
has proven to work well for many projects. When in doubt, a good rule
of thumb is: **Never commit to master**. That means that all 
development should take place within feature branches, and the 
local master branch is updated only by pulling from the upstream
repository.


Core Modules
=============

Most developers will be primarily interested in developing modules 
for the ``shake`` program.
The source for these modules may be found in the directory
*shakemap/coremods*. All of the core modules consist of classes that
inherit from the CoreModule abstract class found in *base.py*. When
developing a new module class, the 
developer must set the class variable ``command_name``, and the 
docstring for the new class should specify this command name followed 
by a brief, one-line description of the module's function. The 
developer must then define the ``execute`` function to perform the 
action of the new module. The event ID will be found in
``self._eventid``. The docstring for the ``execute`` module should
be a more substantial explanation of the module's function and outputs
than is found in the class docstring. See the source for *contour.py*
or *info.py* for examples of the way core modules are implemented.

If a module is properly implemented, the ``shake`` program will discover
it automatically and include it in the list of available modules.


Logging
========
Modules should log to ``self.logger``. We encourage logging generally
useful information as ``self.logger.info``, particularly anything that
might normally be put in as a ``print`` statement. We also encourage
logging copious amounts of potentially useful information as
``logger.debug``. Examples include from the ``model`` core module
include:

.. code-block:: python

    self.logger.info('Inside model')
    self.logger.info('%s: nom bias %f nom stddev %f; %d stations (time=%f sec)'
                     % (imtstr, nominal_bias[imtstr], np.sqrt(nom_variance),
                        np.size(sta_lons_rad[imtstr]), bias_time))
    self.logger.debug('\ttime for %s distance=%f' % (imtstr, ddtime))
    self.logger.debug('\ttime for %s correlation=%f' % (imtstr, ctime))
	


