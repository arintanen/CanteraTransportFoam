/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2024 OpenFOAM Foundation
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
    Foam::functionObjects::specieFlux

Description
    These functions calculate the specie-flux and write it as a
    surfaceScalarField called 'specie<Type>Flux(<specieName>)'. There are
    three such functions; specieAdvectiveFlux and specieDiffusiveFlux return
    the advective and diffusive parts of the flux, respectively, and
    specieFlux returns the total combined flux.

    Example of function object specification:
    \verbatim
    specieFlux
    {
        type    specieFlux; // specieAdvectiveFlux, specieDiffusiveFlux
        libs    ("libfieldFunctionObjects.so");
        field   NH3;
    }
    \endverbatim

    Or, using the standard configuration:
    \verbatim
    #includeFunc specieFlux(NH3)
    \endverbatim

Usage
    \table
        Property     | Description                | Required | Default value
        type         | Type name: specieFlux, \\
                       specieAdvectiveFlux, or \\
                       specieDiffusiveFlux        | yes      |
        field        | Name of the specie/mass \\
                       fraction field             | yes      |
        schemesField | Name of the field from \\
                       which schemes are taken    | no       | Yi
        region       | Region to be evaluated     | no       | default region
    \endtable

SourceFiles
    specieFlux.C

\*---------------------------------------------------------------------------*/

#ifndef specieDiv_H
#define specieDiv_H

#include "fieldExpression.H"
#include "dictionary.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace functionObjects
{


/*---------------------------------------------------------------------------*\
                       Class specieDiffusiveFlux Declaration
\*---------------------------------------------------------------------------*/

class specieDiv
:
    public fieldExpression
{
private:

    // Private Member Functions

        //- Calculate the specie flux field
 	virtual bool calc();

        dictionary specieFluxDict;
        
public:

    //- Runtime type information
    TypeName("specieDiv");


    // Constructors

        //- Construct from Time and dictionary
        specieDiv
        (
            const word& name,
            const Time& runTime,
            const dictionary& dict
        );
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace functionObjects
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
