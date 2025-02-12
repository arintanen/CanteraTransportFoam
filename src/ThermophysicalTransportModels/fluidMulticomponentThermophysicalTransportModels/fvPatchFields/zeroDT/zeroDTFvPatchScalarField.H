/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2023 OpenFOAM Foundation
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
    Foam::zeroDTFvPatchScalarField

Description


Usage


See also
    Foam::fixedValueFvPatchVectorField

SourceFiles
    zeroDTFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef zeroDTFvPatchScalarField_H
#define zeroDTFvPatchScalarField_H

#include "zeroGradientFvPatchFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class zeroDTFvPatchScalarField Declaration
\*---------------------------------------------------------------------------*/

class zeroDTFvPatchScalarField
:
    public zeroGradientFvPatchScalarField
{

public:

    //- Runtime type information
    TypeName("zeroDT");


    // Constructors

        //- Construct from patch and internal field
        zeroDTFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&
        );

        //- Construct from patch, internal field and dictionary
        zeroDTFvPatchScalarField
        (
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const dictionary&
        );

        //- Construct by mapping given zeroDTFvPatchScalarField
        //  onto a new patch
        zeroDTFvPatchScalarField
        (
            const zeroDTFvPatchScalarField&,
            const fvPatch&,
            const DimensionedField<scalar, volMesh>&,
            const fieldMapper&
        );
        //- Copy constructor setting internal field reference
        zeroDTFvPatchScalarField
        (
            const zeroDTFvPatchScalarField&,
            const DimensionedField<scalar, volMesh>&
        );
        
        //- Disallow copy without setting internal field reference
        zeroDTFvPatchScalarField
        (
            const zeroDTFvPatchScalarField&
        ) = delete;



        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone
        (
            const DimensionedField<scalar, volMesh>& iF
        ) const
        {
            return tmp<fvPatchScalarField>
            (
                new zeroDTFvPatchScalarField(*this, iF)
            );
        }


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
