/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "oscillatingSwirlInletVelocityFvPatchVectorField.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "surfaceFields.H"
#include "mathematicalConstants.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

Foam::scalar Foam::oscillatingSwirlInletVelocityFvPatchVectorField::t() const
{
    return db().time().timeOutputValue();
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::oscillatingSwirlInletVelocityFvPatchVectorField::
oscillatingSwirlInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(p, iF),
    rpm_(),
    freq_(),
    Umax_(),
    cavityRadius_(Zero),
    axis_(Zero),
    y_(Zero)
{}


Foam::oscillatingSwirlInletVelocityFvPatchVectorField::
oscillatingSwirlInletVelocityFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchVectorField(p, iF),
    rpm_(readScalar(dict.lookup("rpm"))),
    freq_(readScalar(dict.lookup("freq"))),
    Umax_(readScalar(dict.lookup("Umax"))),
    cavityRadius_(pTraits<vector>(dict.lookup("cavityRadius"))),
    axis_(pTraits<vector>(dict.lookup("axis"))),
    y_(pTraits<vector>(dict.lookup("y")))
{}


Foam::oscillatingSwirlInletVelocityFvPatchVectorField::
oscillatingSwirlInletVelocityFvPatchVectorField
(
    const oscillatingSwirlInletVelocityFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchVectorField(ptf, p, iF, mapper),
    rpm_(ptf.rpm_),
    freq_(ptf.freq_),
    Umax_(ptf.Umax_),
    cavityRadius_(ptf.cavityRadius_),
    axis_(ptf.axis_),
    y_(ptf.y_)
{}


Foam::oscillatingSwirlInletVelocityFvPatchVectorField::
oscillatingSwirlInletVelocityFvPatchVectorField
(
    const oscillatingSwirlInletVelocityFvPatchVectorField& ptf
)
:
    fixedValueFvPatchVectorField(ptf),
    rpm_(ptf.rpm_),
    freq_(ptf.freq_),
    Umax_(ptf.Umax_),
    cavityRadius_(ptf.cavityRadius_),
    axis_(ptf.axis_),
    y_(ptf.y_)
{}


Foam::oscillatingSwirlInletVelocityFvPatchVectorField::
oscillatingSwirlInletVelocityFvPatchVectorField
(
    const oscillatingSwirlInletVelocityFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchVectorField(ptf, iF),
    rpm_(ptf.rpm_),
    freq_(ptf.freq_),
    Umax_(ptf.Umax_),
    cavityRadius_(ptf.cavityRadius_),
    axis_(ptf.axis_),
    y_(ptf.y_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::oscillatingSwirlInletVelocityFvPatchVectorField::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchVectorField::autoMap(m);

}


void Foam::oscillatingSwirlInletVelocityFvPatchVectorField::rmap
(
    const fvPatchVectorField& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchVectorField::rmap(ptf, addr);

    const oscillatingSwirlInletVelocityFvPatchVectorField& tiptf =
        refCast<const oscillatingSwirlInletVelocityFvPatchVectorField>(ptf);


}


void Foam::oscillatingSwirlInletVelocityFvPatchVectorField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    const scalar totArea = gSum(patch().magSf());

    const vector avgCenter = gSum(patch().Cf()*patch().magSf())/totArea;
    const vector avgNormal = gSum(patch().Sf())/totArea;

    // const scalar Umax = Umax_->value(t);
    const scalar sineval = std::sin(2*(constant::mathematical::pi)*freq_*t());

    const vectorField& c = patch().Cf();

    scalarField coord = (c & y_)/(cavityRadius_ & y_);


    // Update angular velocity - convert [rpm] to [rad/s]
    tmp<vectorField> tangentialVelocity
        (
            (rpm_*constant::mathematical::pi/30.0)
          * (patch().Cf() - avgCenter) ^ avgNormal
        );

    tmp<vectorField> n = patch().nf();

    if (sineval>=0)
    {
        // swirl + pulsations during expulsion phase
        fixedValueFvPatchVectorField::operator==
        (
          tangentialVelocity + axis_*Umax_*sineval*(1.0 - sqr(coord))
        );
    }
    else
    {
        // no swirl during suction phase
         fixedValueFvPatchVectorField::operator==
         (
           (axis_*Umax_*sineval)*(1-sqr(coord))
         );
    }

    fixedValueFvPatchField<vector>::updateCoeffs();
}


void Foam::oscillatingSwirlInletVelocityFvPatchVectorField::write
(
    Ostream& os
) const
{
    fvPatchVectorField::write(os);
    os.writeKeyword("rpm") << rpm_ << token::END_STATEMENT << nl;
    os.writeKeyword("freq") << freq_ << token::END_STATEMENT << nl;
    os.writeKeyword("Umax") << Umax_ << token::END_STATEMENT << nl;
    os.writeKeyword("cavityRadius") << cavityRadius_ << token::END_STATEMENT <<nl;
    os.writeKeyword("axis") << axis_ << token::END_STATEMENT << nl;
    writeEntry("value", os);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
   makePatchTypeField
   (
       fvPatchVectorField,
       oscillatingSwirlInletVelocityFvPatchVectorField
   );
}


// ************************************************************************* //
