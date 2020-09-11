/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2012 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::oscillatingSwirlInletVelocityFvPatchVectorField

Group
    grpInletBoundaryConditions

Description
    This boundary condition provides a volumetric- OR mass-flow normal vector
    boundary condition by its magnitude as an integral over its area with a
    swirl component determined by the angular speed, given in revolutions per
    minute (RPM)

    The basis of the patch (volumetric or mass) is determined by the
    dimensions of the flux, phi. The current density is used to correct the
    velocity when applying the mass basis.

    \heading Patch usage

    \table
        Property     | Description             | Required    | Default value
        phi          | flux field name         | no          | phi
        rho          | density field name      | no          | rho
        flowRate     | flow rate profile       | yes         |
        rpm          | rotational speed profile | yes        |
    \endtable

    Example of the boundary condition specification:
    \verbatim
    myPatch
    {
        type             oscillatingSwirlInletVelocity;
        freq             10;
        Umax             0.2;
        rpm              100;
        cavityRadius     (22.5 0 0);
        axis             (0 1 0);
        y                (1 0 0);
        values           (0 0 0);
    }
    \endverbatim

Note
    - the \c flowRate and \c rpm entries are DataEntry types, able to describe
      time varying functions.  The example above gives the usage for supplying
      constant values.
    - the value is positive into the domain

SeeAlso
    Foam::fixedValueFvPatchField

SourceFiles
    oscillatingSwirlInletVelocityFvPatchVectorField.C

\*---------------------------------------------------------------------------*/

#ifndef oscillatingSwirlInletVelocityFvPatchVectorField_H
#define oscillatingSwirlInletVelocityFvPatchVectorField_H

#include "fixedValueFvPatchFields.H"
#include "Function1.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
   Class oscillatingSwirlInletVelocityFvPatchVectorField Declaration
\*---------------------------------------------------------------------------*/

class oscillatingSwirlInletVelocityFvPatchVectorField
:
    public fixedValueFvPatchVectorField
{
    // Private data

        //- Angular speed in revolutions per minute (RPM)
        scalar rpm_;

        //- Frequency of oscillation of the diaphragm
        scalar freq_;

        //- Fully developed flow condition max velocity
        scalar Umax_;

        //- Inlet integral flow rate
        vector cavityRadius_;

        //- Normal vector to the patch (user specified)
        vector axis_;

        //- Transverse Direction
        vector y_;



      //- Private Member Functions

          //- Return current time
          scalar t() const;


public:

   //- Runtime type information
   TypeName("oscillatingSwirlInletVelocity");


   // Constructors

        //- Construct from patch and internal field
        oscillatingSwirlInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        oscillatingSwirlInletVelocityFvPatchVectorField
        (
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given
        //  oscillatingSwirlInletVelocityFvPatchVectorField
        //  onto a new patch
        oscillatingSwirlInletVelocityFvPatchVectorField
        (
            const oscillatingSwirlInletVelocityFvPatchVectorField&,
            const fvPatch&,
            const DimensionedField<vector, volMesh>&,
            const fvPatchFieldMapper&
        );

        //- Construct as copy
        oscillatingSwirlInletVelocityFvPatchVectorField
        (
            const oscillatingSwirlInletVelocityFvPatchVectorField&
        );

        //- Construct and return a clone
        virtual tmp<fvPatchVectorField> clone() const
        {
            return tmp<fvPatchVectorField>
            (
                new oscillatingSwirlInletVelocityFvPatchVectorField(*this)
            );
        }

        //- Construct as copy setting internal field reference
        oscillatingSwirlInletVelocityFvPatchVectorField
        (
            const oscillatingSwirlInletVelocityFvPatchVectorField&,
            const DimensionedField<vector, volMesh>&
        );

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchVectorField> clone
        (
            const DimensionedField<vector, volMesh>& iF
        ) const
        {
            return tmp<fvPatchVectorField>
            (
                new oscillatingSwirlInletVelocityFvPatchVectorField(*this, iF)
            );
        }


    // Member functions

      // Mapping functions

        //- Map (and resize as needed) from self given a mapping object
        virtual void autoMap
        (
            const fvPatchFieldMapper&
        );

        //- Reverse map the given fvPatchField onto this fvPatchField
        virtual void rmap
        (
            const fvPatchVectorField&,
            const labelList&
        );

        // Evaluation

            //- Update the coefficients associated with the patch field
            virtual void updateCoeffs();

            //- Write
        virtual void write(Ostream&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
