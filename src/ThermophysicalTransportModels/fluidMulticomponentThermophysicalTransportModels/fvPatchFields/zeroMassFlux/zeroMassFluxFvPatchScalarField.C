/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2024 OpenFOAM Foundation
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

#include "zeroMassFluxFvPatchScalarField.H"
#include "fieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::zeroMassFluxFvPatchScalarField::zeroMassFluxFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedGradientFvPatchScalarField(p, iF, dict, false),
    curTimeIndex_(-1)
{

    fvPatchField<scalar>::operator=(patchInternalField());
    gradient() = Zero;
    
}


Foam::zeroMassFluxFvPatchScalarField::zeroMassFluxFvPatchScalarField
(
    const zeroMassFluxFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedGradientFvPatchScalarField(ptf, p, iF, mapper, false),
    curTimeIndex_(-1)
{
    map(refCast<const zeroMassFluxFvPatchScalarField>(ptf), mapper);
}


Foam::zeroMassFluxFvPatchScalarField::zeroMassFluxFvPatchScalarField
(
    const zeroMassFluxFvPatchScalarField& wbppsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedGradientFvPatchScalarField(wbppsf, iF),
    curTimeIndex_(-1)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::zeroMassFluxFvPatchScalarField::map
(
    const zeroMassFluxFvPatchScalarField& ptf,
    const fieldMapper& mapper
)
{
    mapper(*this, ptf, [&](){ return patchInternalField(); });
    mapper(gradient(), ptf.gradient(), scalar(0));
}


void Foam::zeroMassFluxFvPatchScalarField::map
(
    const fvPatchScalarField& ptf,
    const fieldMapper& mapper
)
{
    map(refCast<const zeroMassFluxFvPatchScalarField>(ptf), mapper);
}


void Foam::zeroMassFluxFvPatchScalarField::updateCoeffs
(
    const scalarField& snGradp
)
{
    if (updated())
    {
        return;
    }
    curTimeIndex_ = this->db().time().timeIndex();

    gradient() = snGradp;
    fixedGradientFvPatchScalarField::updateCoeffs();
}


void Foam::zeroMassFluxFvPatchScalarField::updateCoeffs()
{
    if (updated())
    {
        return;
    }

    if (curTimeIndex_ != this->db().time().timeIndex())
    {
        FatalErrorInFunction
            << "updateCoeffs(const scalarField& snGradp) MUST be called before"
               " updateCoeffs() or evaluate() to set the boundary gradient."
            << exit(FatalError);
    }
}


void Foam::zeroMassFluxFvPatchScalarField::write(Ostream& os) const
{
    fixedGradientFvPatchScalarField::write(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        zeroMassFluxFvPatchScalarField
    );
}


// ************************************************************************* //
