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

#include "specieFlux2.H"
#include "fluidThermophysicalTransportModel.H"
#include "fvcFlux.H"
#include "fvcSnGrad.H"
#include "multicomponentThermo.H"
#include "addToRunTimeSelectionTable.H"



// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace functionObjects
{
    defineTypeNameAndDebug(specieSoretFlux, 0);
    addToRunTimeSelectionTable(functionObject, specieSoretFlux, dictionary);

    defineTypeNameAndDebug(specieCorrectionFlux, 0);
    addToRunTimeSelectionTable(functionObject, specieCorrectionFlux, dictionary);

    defineTypeNameAndDebug(specieOrdinaryDiffusionFlux, 0);
    addToRunTimeSelectionTable(functionObject, specieOrdinaryDiffusionFlux, dictionary);
}
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //



Foam::tmp<Foam::surfaceScalarField>
Foam::functionObjects::specieSoretFlux::calc
(
    const fluidThermophysicalTransportModel& ttm,
    const volScalarField& Yi
)
{

    const volScalarField & DTi = mesh_.lookupObject<volScalarField>("DT_" + Yi.name());
    const volScalarField & T = ttm.thermo().T();
    return -fvc::interpolate(DTi)*fvc::snGrad(T)/fvc::interpolate(T);
}


Foam::tmp<Foam::surfaceScalarField>
Foam::functionObjects::specieCorrectionFlux::calc
(
    const fluidThermophysicalTransportModel& ttm,
    const volScalarField& Yi
)
{
    const surfaceScalarField & Jc = mesh_.lookupObject<surfaceScalarField>("Jc");
    return -Jc*fvc::interpolate(Yi);
}


Foam::tmp<Foam::surfaceScalarField>
Foam::functionObjects::specieOrdinaryDiffusionFlux::calc
(
    const fluidThermophysicalTransportModel& ttm,
    const volScalarField& Yi
)
{
    const volScalarField & DTi = mesh_.lookupObject<volScalarField>("DT_" + Yi.name());
    const volScalarField & T = ttm.thermo().T();
    return ttm.j(Yi) + fvc::interpolate(DTi)*fvc::snGrad(T)/fvc::interpolate(T);
}


// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //



// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //



Foam::functionObjects::specieSoretFlux::specieSoretFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    specieFluxBase(name, runTime, dict, typeName)
{}


Foam::functionObjects::specieCorrectionFlux::specieCorrectionFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    specieFluxBase(name, runTime, dict, typeName)
{}


Foam::functionObjects::specieOrdinaryDiffusionFlux::specieOrdinaryDiffusionFlux
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    specieFluxBase(name, runTime, dict, typeName)
{}



// ************************************************************************* //
