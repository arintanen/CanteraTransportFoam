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
    Foam::turbulenceThermophysicalTransportModels::CanteraMulticomponent

Description
    Base class for Fickian transport
    
SourceFiles
    CanteraMulticomponent.C

\*---------------------------------------------------------------------------*/


#ifndef CanteraMulticomponent_H
#define CanteraMulticomponent_H
#include "cantera/core.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#include "LUscalarMatrix.H"
#include "Function2.H"
namespace Foam
{
namespace laminarThermophysicalTransportModels
{

/*---------------------------------------------------------------------------*\
                Class CanteraMulticomponent Declaration
\*---------------------------------------------------------------------------*/

template<class laminarThermophysicalTransportModel>
class CanteraMulticomponent
:
    public unityLewisFourier<laminarThermophysicalTransportModel>,
    public TopoChangeableMeshObject<fvMesh>
{
    // Private data
            
	   	    
        //- Enable soret effect  
	bool soretEnabled_; 

            //- Generalised Fick's law diffusion coefficients field list.
            //  This is the diagonal of the mass diffusion coefficient matrix
            mutable PtrList<volScalarField> Dii_;

            //- List of fields of the explicit part of the mass diffusion flux
            // of the species
            mutable PtrList<surfaceScalarField> jexp_;


        // Workspace for diffusion coefficient transformation

            //- Molecular weights of the species
            mutable scalarField W;

            //- List of mass-fraction field pointers
            //  for the internal mesh or a patch
            mutable List<const scalarField*> YPtrs;

            //- Matrix of mass diffusion coefficient field pointers
            //  for the internal mesh or a patch
            mutable SquareMatrix<scalarField*> DijPtrs;

            //- Mass-fractions at a cell or face
            mutable scalarField Y;

            //- Mole-fractions at a cell or face
            mutable scalarField X;

            //- Binary mass diffusion coefficient matrix at a cell or face
            mutable scalarSquareMatrix DD;

            //- Matrix form of the coefficients in the Maxwell-Stefan equation
            //  at a cell or face
            mutable LUscalarMatrix A;

            //- Matrix of composition coefficients at a cell or face
            //  used to transform the binary mass diffusion coefficients into
            //  the generalised Fick's law diffusion coefficients
            mutable scalarSquareMatrix B;

            //- Inverse of A
            mutable scalarSquareMatrix invA;

            //- Matrix of the generalised Fick's law diffusion coefficients
            //  at a cell or face
            mutable scalarSquareMatrix D;


    // Private member functions

        //- Transform the Binary mass diffusion coefficient matrix DD into the
        //  matrix of the generalised Fick's law diffusion coefficients D
        void transformDiffusionCoefficient() const;

        //- Calls transformDiffusionCoefficient() for each cell and patch face
        //  and maps the resulting D into DijPtrs
        void transformDiffusionCoefficientFields() const;

        //- Calls transformDiffusionCoefficientFields() for the internal mesh
        //  and each patch and maps the resulting DijPtrs into Dii_ and Dij
        void transform(List<PtrList<volScalarField>>& Dij) const;

        //- Update Dii_
        void updateDii() const;

        //- Update Dii_ if not yet constructed and return
        const PtrList<volScalarField>& Dii() const;

        //- Update jexp_ if not yet constructed and return
        const PtrList<surfaceScalarField>& jexp() const;

        //- Thermal diffusion coefficients in the mixture 
        mutable PtrList<volScalarField> DT_;    
            //- Array of specie binary mass diffusion coefficient functions
            //  [m^2/s]
            List<PtrList<Function2<scalar>>> DFuncs_;        
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
        CanteraMulticomponent
        (
            const momentumTransportModel& momentumTransport,
            const thermoModel& thermo
        );

        TypeName("CanteraMulticomponent");
    //- Destructor
    virtual ~CanteraMulticomponent()
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
    #include "CanteraMulticomponent.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
