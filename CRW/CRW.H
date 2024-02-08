/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2020 OpenFOAM Foundation
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
    Foam::CRW

Description
    The velocity is perturbed in random direction, with a
    Gaussian random number distribution with variance sigma.
    where sigma is defined below

\*---------------------------------------------------------------------------*/

#ifndef CRW_H
#define CRW_H

#include "DispersionRASModel.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                   Class CRW Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CRW
:
    public DispersionRASModel<CloudType>
{
protected:

    // Protected data

        // Locally cached turbulence fields

            //- Gradient of sigma
            const volVectorField* gradSigmaPtr_;

            //- Take ownership of the grad(sigma)
            mutable bool ownGradSigma_;
            
            //- Gradient of sigmaCL
            const volVectorField* gradSigmaCLPtr_;

            //- Take ownership of the grad(sigmaCL)
            mutable bool ownGradSigmaCL_;
            
            //- Gradient of sigmaCL
            const volVectorField* gradSigmaCL4Ptr_;

            //- Take ownership of the grad(sigmaCL)
            mutable bool ownGradSigmaCL4_;
            
            
            
public:

    //- Runtime type information
    TypeName("CRW");

    
    // Constructors

        //- Construct from components
        CRW(const dictionary& dict, CloudType& owner);

        //- Construct copy
        CRW(const CRW<CloudType>& dm);

        //- Construct and return a clone
        virtual autoPtr<DispersionModel<CloudType>> clone() const
        {
            return autoPtr<DispersionModel<CloudType>>
            (
                new CRW<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~CRW();


    // Member Functions

	//- Cache carrier fields
        virtual void cacheFields(const bool store);

        //- Update (disperse particles)
        virtual vector update
        (
            //const typename CloudType::parcelType& p,
               const scalar dt,
    		const label celli,
    		const vector& U,
    		const vector& Uc,
    		vector& UTurb,
    		vector& UTurbn,
    		vector& RMS,
    		vector& RMSn,
    		scalar& tTurb
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CRW.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
