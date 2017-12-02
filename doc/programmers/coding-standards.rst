Coding standards
================

General Object-Oriented design principles you should try to follow:
  1. Identify the aspects of your application that vary and separate them from what stays the same;
  2. Program to an interface, not an implementation;
  3. Favor composition over inheritance;
  4. Strive for loosely coupled designs between objects that interact;
  5. Classes should be open for extension, but closed for modification;
  6. Depend upon abstractions. Do not depend upon concrete classes;
  7. Principle of Least Knowledge. Talk only to your immediate friends;

:cite:`Sutter2004,Cline1998,CppFAQs`

Including header files
----------------------

Do not include header files unnecessarily. Unnecessary include
directives and/or forward declarations introduce nasty
interdependencies among different parts of the code.  This reflects
mainly in longer compilation times, but also in uglier looking code
(see also the discussion in :cite:`Sutter1999`).

Follow these guidelines to decide whether to include or forward declare:
  1. class A makes no reference to class B. Neither include nor forward declare B;
  2. class A refers to class B as a friend. Neither include nor forward declare B;
  3. class A contains a pointer/reference to a class B object. Forward declare B;
  4. class A contains functions with a class B object (value/pointer/reference) as parameter/return value. Forward declare B;
  5. class A is derived from class B. include B;
  6. class A contains a class B object. include B.

.. code-block:: cpp

    #pragma once

    //==============================
    // Forward declared dependencies
    class Foo;
    class Bar;

    //==============================
    // Included dependencies
    #include <vector>
    #include "Parent.hpp"

    //==============================
    // The actual class
    class MyClass : public Parent // Parent object, so #include "Parent.h"
    {
      public:
        std::vector<int> avector; // vector object, so #include <vector>
        Foo * foo;                // Foo pointer, so forward declare
        void Func(Bar & bar);     // Bar reference as parameter, so forward declare

        friend class MyFriend;    // friend declaration is not a dependency
                                  //    don't do anything about MyFriend
    };


Proper overloading of `operator<<`
----------------------------------

Suppose we have an inheritance hierarchy made of an abstract base class, Base, and
two derived classes, Derived1 and Derived2.
In the Base class header file we will define a pure virtual private function printObject
and provide a public friend overload of operator<<:

.. code-block:: cpp

    #include <iosfwd>

    class Base
    {
      public:
        // All your other very fancy public members
        friend std::ostream & operator<<(std::ostream & os, Base & base)
        {
                return base.printObject(os);
        }
      protected:
        // All your other very fancy protected members
      private:
        // All your other very fancy private members
        virtual std::ostream & printObject(std::ostream & os) = 0;
    }

The printObject method can also be made (impure) virtual, it really depends on your class hierarchy.
Derived1 and Derived2 header files will provide a public friend overload of operator<< (friendliness
isn't inherited, transitive or reciprocal) and an override for the printObject method:

.. code-block:: cpp

    #include <iosfwd>

    #include "Base.hpp"

    class Derived1 : public Base
    {
      public:
        // All your other very fancy public members
        friend std::ostream & operator<<(std::ostream & os, Derived1 & derived)
        {
          return derived.printObject(os);
        }
      protected:
        // All your other very fancy protected members
      private:
        // All your other very fancy private members
        virtual std::ostream & printObject(std::ostream & os);
    }

    class Derived2 : public Base
    {
      public:
        // All your other very fancy public members
        friend std::ostream & operator<<(std::ostream & os, Derived2 & derived)
        {
          return derived.printObject(os);
        }
      protected:
        // All your other very fancy protected members
      private:
        // All your other very fancy private members
        virtual std::ostream & printObject(std::ostream & os);
    }

Code formatting
---------------

We conform to the so-called Linux (aka kernel) formatting style for C/C++ code
(see http://en.wikipedia.org/wiki/Indent_style#Kernel_style) with minimal
modifications.
Using `clang-format <https://clang.llvm.org/docs/ClangFormat.html>`_ is the
preferred method to get the source code in the right format.
Formatting style is defined in the ``.clang-format`` file, kept at the root of the project.

.. note::
   We recommend using at least v3.9 of the program, which is the version used to
   generate the ``.clang-format`` file defining all formatting settings.

``clang-format`` can be `integrated with both
Emacs and Vim. <https://clang.llvm.org/docs/ClangFormat.html#vim-integration>`_
It is also possible to install the Git pre-commit hooks to perform the necessary code style
checks prior to committing changes:

.. code-block:: bash

   cd .git/hooks
   cp --symbolic-link ../../.githooks/* .
