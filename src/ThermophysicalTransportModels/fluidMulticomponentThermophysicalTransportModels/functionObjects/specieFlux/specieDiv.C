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

\*---------------------------------------------------------------------------*/

#include "specieDiv.H"
#include "specieFlux.H"
#include "fluidThermophysicalTransportModel.H"
#include "fvcFlux.H"
#include "fvcDiv.H"
#include "fvcSnGrad.H"
#include "multicomponentThermo.H"
#include "addToRunTimeSelectionTable.H"



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(specieDiv, 0);
    addToRunTimeSelectionTable(functionObject, specieDiv, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //





// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

bool Foam::functionObjects::specieDiv::calc()
{

    specieFlux sFlux("specieFlux", mesh_.time(), specieFluxDict);
    sFlux.execute();
    const surfaceScalarField& fluxField =
        lookupObject<surfaceScalarField>(resultName_ + "flux");
        
    
    store(resultName_, fvc::div(mesh_.magSf()*fluxField));
    sFlux.clear();
    
    return true;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //



Foam::functionObjects::specieDiv::specieDiv
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    fieldExpression(name, runTime, dict, typeName)
{
    word fluxType = dict.lookup("fluxType");
    if(!dict.found("result")) 
    {
        resultName_ = word("div" + fluxType + '(' + fieldName_ + ')');   
    }
    specieFluxDict.add("field", fieldName_);
    specieFluxDict.add("fluxType", fluxType);
    specieFluxDict.add("result", resultName_ + "flux");
}



// ************************************************************************* //
