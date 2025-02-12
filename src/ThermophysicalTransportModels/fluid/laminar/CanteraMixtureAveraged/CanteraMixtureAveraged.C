/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
Copyright (C) 2024 Aleksi Rintanen, Aalto university
License

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

#include "CanteraMixtureAveraged.H"
#include "fvcDiv.H"
#include "fvcLaplacian.H"
#include "fvcSnGrad.H"
#include "fvmSup.H"
#include "fvmDiv.H"
#include "surfaceInterpolate.H"
#include "zeroDTFvPatchScalarField.H"
#include "zeroMassFluxFvPatchScalarField.H"
#include "cantera/transport.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarThermophysicalTransportModels
{
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class laminarThermophysicalTransportModel>
CanteraMixtureAveraged<laminarThermophysicalTransportModel>::CanteraMixtureAveraged
(
    const momentumTransportModel& momentumTransport,
    const thermoModel& thermo
)
:
    unityLewisFourier<laminarThermophysicalTransportModel>
    (
        typeName,
        momentumTransport,
        thermo
    ),
    Jc_
     (
         IOobject
         (
             thermo.phasePropertyName("Jc"),
             this->thermo().T().mesh().time().name(),
             this->thermo().T().mesh(),
             IOobject::NO_READ,
             IOobject::NO_WRITE
         ),
         this->thermo().T().mesh(),
         dimensionedScalar(dimMass/dimArea/dimTime, 0)
     ),

    soretEnabled_(true)
    
{
    if(this->thermo().he().name() == "e") {
        FatalErrorInFunction
            << "Internal energy is not supported as a solution variable. Please use enthalpy instead."
            << exit(FatalError);
    }

    const PtrList<volScalarField>& Y = this->thermo().Y();
    DT_.setSize(Y.size());
    forAll(DT_, i)
     {
         DT_.set
         (
             i,
             new volScalarField
             (
                 IOobject
                 (
                     "DT_" + Y[i].name(),
                     this->thermo().mesh().time().name(),
                     this->thermo().mesh(),
                     IOobject::NO_READ,
                     IOobject::NO_WRITE
                 ),
                 this->thermo().mesh(),
                 dimensionedScalar(dimDynamicViscosity, 0)
             )
         );

     }
    read();
    this->printCoeffs(typeName);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class laminarThermophysicalTransportModel>
bool CanteraMixtureAveraged<laminarThermophysicalTransportModel>::read()
{
    if
    (
        unityLewisFourier<laminarThermophysicalTransportModel>::read()
    )
    {
        soretEnabled_ = this->coeffDict_.lookupOrDefault("soretEnabled",true);
        Info << "Reading Cantera mechanism" << endl;
        fileName mechname = this->coeffDict_.lookup("mechFile");
        solution_ = Cantera::newSolution(mechname);	
        const PtrList<volScalarField>& Y = this->thermo().Y();
        speciesIndex_.setSize(Y.size());

        if(Y.size() != static_cast<label>(solution_->thermo()->nSpecies()))
        {
            FatalErrorInFunction
                << "Openfoam and Cantera mechanism have different number of species."
                << exit(FatalError);                
        }
        // Cantera uses std containers
        const std::vector<std::string> & speciesNames = solution_->thermo()->speciesNames();
        forAll(Y, i) 
        {
            const std::vector<std::string>::const_iterator iter = std::find(speciesNames.begin(), speciesNames.end(), Y[i].name());
            if(iter == speciesNames.end()) 
            {
                FatalErrorInFunction
                    << "Specie " << Y[i].name() << " not found in Cantera mechanism."
                    << exit(FatalError);                
            }
            speciesIndex_[i] = std::distance(speciesNames.begin(), iter);
        }
        updateDm();
        return true;
    }
    else
    {
        return false;
    }
}



template<class laminarThermophysicalTransportModel>
tmp<surfaceScalarField> CanteraMixtureAveraged<laminarThermophysicalTransportModel>::q() const
{
    tmp<surfaceScalarField> tmpq
    (
        surfaceScalarField::New
        (
            IOobject::groupName
            (
                "q",
                this->momentumTransport().alphaRhoPhi().group()
            ),
           -fvc::interpolate(this->alpha()*this->kappaEff())
           *fvc::snGrad(this->thermo().T())
        )
    );

    const PtrList<volScalarField>& Y = this->thermo().Y();
    const volScalarField& p = this->thermo().p();
    const volScalarField& T = this->thermo().T();
    
    if (Y.size())
    {
        surfaceScalarField sumJh
        (
            surfaceScalarField::New
            (
                "sumJh",
                Y[0].mesh(),
                dimensionedScalar(dimMass/dimArea/dimTime*dimEnergy/dimMass, 0)
            )
        );

        forAll(Y, i)
        {

            const volScalarField hi(this->thermo().hsi(i, p, T));
            const surfaceScalarField ji(unityLewisFourier<laminarThermophysicalTransportModel>::j(Y[i]));
            sumJh += ji*fvc::interpolate(hi);
            
        }

        tmpq.ref() += sumJh;
    }
    tmpq.ref() -= Jc_*fvc::interpolate(this->thermo().he());
    return tmpq;
}


template<class laminarThermophysicalTransportModel>
tmp<fvScalarMatrix> CanteraMixtureAveraged<laminarThermophysicalTransportModel>::divq
(
    volScalarField& he
) const
{
     correctJc();
     
     tmp<fvScalarMatrix> tmpDivq
     (
        fvm::Su 
         (
             -fvc::laplacian(this->alpha()*this->kappaEff(), this->thermo().T()),
             he
         )
     );
     
     const PtrList<volScalarField>& Y = this->thermo().Y();
     const volScalarField& p = this->thermo().p();
     const volScalarField& T = this->thermo().T();


     
     tmpDivq.ref() -=
        correction(fvm::laplacian(this->alpha()*this->alphaEff(), he));
     

 
     surfaceScalarField sumJh
     (
        surfaceScalarField::New
        (
            "sumJh",
            he.mesh(),
            dimensionedScalar(dimMass/dimArea/dimTime*he.dimensions(), 0)
        )
     );
 
     forAll(Y, i)
     {

        const volScalarField hi(this->thermo().hsi(i, p, T));

        const surfaceScalarField ji(unityLewisFourier<laminarThermophysicalTransportModel>::j(Y[i]));

        sumJh += ji*fvc::interpolate(hi);
     
     }
 
 
    tmpDivq.ref() += fvc::div(sumJh*he.mesh().magSf());
    tmpDivq.ref() -= fvm::div(Jc_*he.mesh().magSf(),he,"div(phi,Yi_h)");
    return tmpDivq;
 }

template<class laminarThermophysicalTransportModel>
tmp<surfaceScalarField> CanteraMixtureAveraged<laminarThermophysicalTransportModel>::j
(
    const volScalarField& Yi
) const
{
        const volScalarField& T = this->thermo().T();
        return unityLewisFourier<laminarThermophysicalTransportModel>::j(Yi)
           - fvc::interpolate
             (
               DT_[this->thermo().specieIndex(Yi)]
             )
            *fvc::snGrad(T)/fvc::interpolate(T)
             - Jc_*fvc::interpolate(Yi); 
}


template<class laminarThermophysicalTransportModel>
tmp<fvScalarMatrix> CanteraMixtureAveraged<laminarThermophysicalTransportModel>::divj
(
    volScalarField& Yi
) const
{
    const volScalarField& T = this->thermo().T();
    return unityLewisFourier<laminarThermophysicalTransportModel>::divj(Yi)
           - fvc::div
             (
                 fvc::interpolate
                 (
                   DT_[this->thermo().specieIndex(Yi)]
                 )
                *fvc::snGrad(T)/fvc::interpolate(T)
                *T.mesh().magSf()
             )
            - fvm::div(Jc_*T.mesh().magSf(),Yi,"div(phi,Yi_h)");  
}


template<class laminarThermophysicalTransportModel>
void CanteraMixtureAveraged<laminarThermophysicalTransportModel>::updateDm() const
{

    const PtrList<volScalarField>& Y = this->thermo().Y();
    const volScalarField& p = this->thermo().p();
    const volScalarField& T = this->thermo().T();
    const volScalarField Wm(this->thermo().W());
    Dm_.setSize(Y.size());
    forAll(Dm_, i)
    {
        Dm_.set
            (
                i,
                volScalarField::New
                (
                    "Dm_" + Y[i].name(),
                    T.mesh(),
                    dimensionedScalar(dimKinematicViscosity, scalar(0))
                )
            );
         DT_[i] = Zero;
        } 
    Cantera::shared_ptr<Cantera::ThermoPhase> gas = solution_->thermo();
    scalarField Ys(Y.size());
    scalarField Dms(Y.size());
    scalarField DTs(Y.size());
    scalar epsilon = 1-Y.size()*SMALL;
    forAll(T,celli) {
    	forAll(Y,i) {
            if(Y[i][celli] > epsilon)
            {
                Ys[speciesIndex_[i]] = epsilon;
            } 
            else {
                Ys[speciesIndex_[i]] = Y[i][celli] + SMALL;
            }
    	}
    	gas->setState_TPY(T[celli], p[celli], Ys.begin());
        solution_->transport()->getMixDiffCoeffs(Dms.begin());
        forAll(Y,i) {
            Dm_[i][celli] = Dms[speciesIndex_[i]];
        }
        if(soretEnabled_) {    
       	    solution_->transport()->getThermalDiffCoeffs(DTs.begin());
	        forAll(Y,i) {
		        DT_[i][celli] = DTs[speciesIndex_[i]];
	        }
       }
    }

    forAll(T.boundaryField(), patchi)
    {
        forAll(T.boundaryField()[patchi],facei) 
        {
            forAll(Y,i) {
                if(Y[i].boundaryField()[patchi][facei] > epsilon)
                {
                    Ys[speciesIndex_[i]] = epsilon;
                } 
                else {
                    Ys[speciesIndex_[i]] = Y[i].boundaryField()[patchi][facei] + SMALL;
                }
    	    }
            gas->setState_TPY(T.boundaryField()[patchi][facei], p.boundaryField()[patchi][facei], Ys.begin());
            solution_->transport()->getMixDiffCoeffs(Dms.begin());
            forAll(Y,i) {
                Dm_[i].boundaryFieldRef()[patchi][facei] = Dms[speciesIndex_[i]];
                
            }
            if(soretEnabled_) {
                solution_->transport()->getThermalDiffCoeffs(DTs.begin());
                forAll(Y,i) {
                    DT_[i].boundaryFieldRef()[patchi][facei] = DTs[speciesIndex_[i]];
                }
            }
        }
        forAll(Y,i)
        {
            if(isA<zeroDTFvPatchScalarField>(Y[i].boundaryField()[patchi]))
            {
                DT_[i].boundaryFieldRef()[patchi] = Zero;          
            }
        }
    }      
      correctJc();
}

 template<class laminarThermophysicalTransportModel>
 void CanteraMixtureAveraged<laminarThermophysicalTransportModel>::predict()
 {
     unityLewisFourier<laminarThermophysicalTransportModel>::predict();
     updateDm();
 }
  
  
 template<class laminarThermophysicalTransportModel>
 bool CanteraMixtureAveraged<laminarThermophysicalTransportModel>::movePoints()
 {
     return true;
 }
  
  
 template<class laminarThermophysicalTransportModel>
 void CanteraMixtureAveraged<laminarThermophysicalTransportModel>::topoChange
 (
     const polyTopoChangeMap& map
 )
 {
     // Delete the cached Dm, will be re-created in predict
     Dm_.clear();
 }
  
  
template<class laminarThermophysicalTransportModel>
void CanteraMixtureAveraged<laminarThermophysicalTransportModel>::mapMesh
 (
     const polyMeshMap& map
 )
 {
     // Delete the cached Dm, will be re-created in predict
     Dm_.clear();
 }
  
  
template<class laminarThermophysicalTransportModel>
void CanteraMixtureAveraged<laminarThermophysicalTransportModel>::distribute
 (
     const polyDistributionMap& map
 )
 {
     // Delete the cached Dm, will be re-created in predict
     Dm_.clear();
 }

 template<class laminarThermophysicalTransportModel>
 const PtrList<volScalarField>&
 CanteraMixtureAveraged<laminarThermophysicalTransportModel>::Dm() const
 {
     if (!Dm_.size())
     {
         updateDm();
     }
  
     return Dm_;
 }
   
template<class laminarThermophysicalTransportModel>
void CanteraMixtureAveraged<laminarThermophysicalTransportModel>::correctJc() const
{
   
   PtrList<volScalarField>& Y = const_cast<PtrList<volScalarField>&>(this->thermo().Y());
   const volScalarField& T = this->thermo().T();
   surfaceScalarField gradLogT = fvc::snGrad(T)/fvc::interpolate(T);
   Jc_ *= scalar(0);
   forAll(Y, i)
   {
   
       surfaceScalarField aDeff = fvc::interpolate(this->alpha()*this->DEff(Y[i]));
       surfaceScalarField DTi = fvc::interpolate(DT_[i]);
       setSnGrad<zeroMassFluxFvPatchScalarField>
            (
                Y[i].boundaryFieldRef(),
                -gradLogT.boundaryField()
               *DTi.boundaryField()
               /aDeff.boundaryField()  
            );
        Jc_ -= aDeff*fvc::snGrad(Y[i]) + DTi*gradLogT;
   }  

}

template<class laminarThermophysicalTransportModel>
tmp<volScalarField> CanteraMixtureAveraged<laminarThermophysicalTransportModel>::DEff
(
    const volScalarField& Yi
) const
{
    return volScalarField::New
    (
        "DEff",
        this->momentumTransport().rho()
       *this->Dm()[this->thermo().specieIndex(Yi)]
    );
}


template<class laminarThermophysicalTransportModel>
tmp<scalarField> CanteraMixtureAveraged<laminarThermophysicalTransportModel>::DEff
(
    const volScalarField& Yi,
    const label patchi
) const
{  

    return
        this->momentumTransport().rho().boundaryField()[patchi]
       *this->Dm()[this->thermo().specieIndex(Yi)].boundaryField()[patchi];
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace laminarThermophysicalTransportModels
} // End namespace Foam

// ************************************************************************* //
