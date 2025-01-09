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
Class
    Foam::turbulenceThermophysicalTransportModels::CanteraMixtureAveraged

Description
    Base class for Fickian transport
    
SourceFiles
    CanteraMixtureAveraged.C

\*---------------------------------------------------------------------------*/


#ifndef CanteraMixtureAveraged_H
#define CanteraMixtureAveraged_H
#include "cantera/core.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace laminarThermophysicalTransportModels
{

/*---------------------------------------------------------------------------*\
                Class CanteraMixtureAveraged Declaration
\*---------------------------------------------------------------------------*/

template<class laminarThermophysicalTransportModel>
class CanteraMixtureAveraged
:
    public unityLewisFourier<laminarThermophysicalTransportModel>
{
    // Private data
            
        //- Correction flux
        mutable surfaceScalarField Jc_;
	   	    
        //- Enable soret effect  
	bool soretEnabled_; 

        //- Calculate correction flux    	    
 	    void correctJc() const;
 	    
        //- Update Dm_
        void updateDm() const;

        //- Mass diffusion coefficients in the mixture
        mutable PtrList<volScalarField> Dm_;
        //- Thermal diffusion coefficients in the mixture 
        mutable PtrList<volScalarField> DT_;    
        
        //- Cantera solution object
        Cantera::shared_ptr<Cantera::Solution> solution_;
        //- Indexing array (if the species are in different order in cantera)
        List<label> speciesIndex_; 	          
protected:

        //- Update Dm_ if not yet constructed and return
        const PtrList<volScalarField>& Dm() const;                
             	    
public:

    typedef typename laminarThermophysicalTransportModel::alphaField
        alphaField;

    typedef typename
        laminarThermophysicalTransportModel::momentumTransportModel
        momentumTransportModel;

    typedef typename laminarThermophysicalTransportModel::thermoModel
        thermoModel;
    // Constructors

        //- Construct from a momentum transport model and a thermo model
        CanteraMixtureAveraged
        (
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo
        );

        TypeName("CanteraMixtureAveraged");
    //- Destructor
    virtual ~CanteraMixtureAveraged()
    {}


    // Member Functions

        //- Read thermophysicalTransport dictionary
        virtual bool read();

        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction [kg/m/s]
        virtual tmp<volScalarField> DEff(const volScalarField& Yi) const;

        //- Effective mass diffusion coefficient
        //  for a given specie mass-fraction for patch [kg/m/s]
        virtual tmp<scalarField> DEff
        (
            const volScalarField& Yi,
            const label patchi
        ) const;
        
        //- Return the heat flux [W/m^2]
        virtual tmp<surfaceScalarField> q() const;

        //- Return the source term for the energy equation
        virtual tmp<fvScalarMatrix> divq(volScalarField& he) const;

        //- Return the specie flux for the given specie mass-fraction [kg/m^2/s]
        virtual tmp<surfaceScalarField> j(const volScalarField& Yi) const;

        //- Return the source term for the given specie mass-fraction equation
        virtual tmp<fvScalarMatrix> divj(volScalarField& Yi) const;

         //- Update the diffusion coefficients
         virtual void predict();
         
           // Mesh changes
  
             //- Update for mesh motion
             virtual bool movePoints();
  
             //- Update topology using the given map
             virtual void topoChange(const polyTopoChangeMap& map);
  
             //- Update from another mesh using the given map
             virtual void mapMesh(const polyMeshMap& map);
  
             //- Redistribute or update using the given distribution map
             virtual void distribute(const polyDistributionMap& map);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
} // End namespace laminarThermophysicalTransportModels
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CanteraMixtureAveraged.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
