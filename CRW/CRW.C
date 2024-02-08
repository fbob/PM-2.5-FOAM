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

\*---------------------------------------------------------------------------*/

#include "CRW.H"
#include "constants.H"
#include "momentumTransportModel.H"
#include "wallDist.H"
#include "demandDrivenData.H"
#include "fvcGrad.H"


using namespace Foam::constant::mathematical;

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CRW<CloudType>::CRW
(
    const dictionary& dict,
    CloudType& owner
)
:
    DispersionRASModel<CloudType>(dict, owner),

    gradSigmaPtr_(nullptr),
    ownGradSigma_(false),
    gradSigmaCLPtr_(nullptr),
    ownGradSigmaCL_(false),
    gradSigmaCL4Ptr_(nullptr),
    ownGradSigmaCL4_(false)
       
{}


template<class CloudType>
Foam::CRW<CloudType>::CRW
(
    const CRW<CloudType>& dm
)
:
    DispersionRASModel<CloudType>(dm),
    gradSigmaPtr_(dm.gradSigmaPtr_),
    ownGradSigma_(dm.ownGradSigma_),
    gradSigmaCLPtr_(dm.gradSigmaCLPtr_),
    ownGradSigmaCL_(dm.ownGradSigmaCL_),
    gradSigmaCL4Ptr_(dm.gradSigmaCL4Ptr_),
    ownGradSigmaCL4_(dm.ownGradSigmaCL4_)
   
{
    dm.ownGradSigma_ = false;
    dm.ownGradSigmaCL_ = false;
    dm.ownGradSigmaCL4_ = false;

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CRW<CloudType>::~CRW()
{
	cacheFields(false);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::CRW<CloudType>::cacheFields(const bool store)
{
    DispersionRASModel<CloudType>::cacheFields(store);

    const objectRegistry& obr = this->owner().mesh();
    const word turbName = IOobject::groupName
    (
         momentumTransportModel::typeName,
         this->owner().U().group()
    );

    
    const volSymmTensorField R_ = obr.lookupObject<volSymmTensorField>("R"); //on récupère le champs tensoriel du tenseur de Reynolds
    
    const volScalarField sigmaXY = R_.component(1);   //on récupère la composante croisée /u'v'
    const volScalarField sigmaY = R_.component(3);    //on récupère la composante /v'² dans la direction normale à la paroi
    
    
/////////////////////Calcul du gradient de k////////////////////////////////////////////////////////////////////////////////// 
    
    if (store)
    {
        gradSigmaPtr_ = fvc::grad(*this->kPtr_).ptr();
        ownGradSigma_ = true;
    }
    else
    {
        if (ownGradSigma_)
        {
            deleteDemandDrivenData(gradSigmaPtr_);
            gradSigmaPtr_ = nullptr;
            ownGradSigma_ = false;
        }
    }
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
/////////////////////Calcul du gradient de la correlation croisée /uv du tenseur de Reynolds/////////////////////////////////////
    
    if (store)
    {
      	gradSigmaCLPtr_ = fvc::grad(sigmaXY/sqrt(sigmaY)).ptr();
      	ownGradSigmaCL_ = true;
    }
    else
    {
        if (ownGradSigmaCL_)
        {
            deleteDemandDrivenData(gradSigmaCLPtr_);
            gradSigmaCLPtr_ = nullptr;
            ownGradSigmaCL_ = false;
        }
    }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
///////////////////Calcul du gradient de la tension de Reynolds dans la direction normale à la paroi///////////////////////////
    
    if (store)
    {
      	gradSigmaCL4Ptr_ = fvc::grad(sqrt(sigmaY)).ptr();
      	ownGradSigmaCL4_ = true;
    }
    else
    {
        if (ownGradSigmaCL4_)
        {
            deleteDemandDrivenData(gradSigmaCL4Ptr_);
            gradSigmaCL4Ptr_ = nullptr;
            ownGradSigmaCL4_ = false;
        }
    }
////////////////////////////////////////////////////////////////////////////////////////////////////////////    
}



template<class CloudType>
Foam::vector Foam::CRW<CloudType>::update
(
    const scalar dt,
    const label celli,
    const vector& U,
    const vector& Uc,
    vector& UTurb,
    vector& UTurbn,
    vector& RMS,
    vector& RMSn,
    scalar& tTurb
)
{
    
///////////Définition de grandeurs nécéssaires pour le calcul ///////////////////////////////////////   
////////// Pour une application plus généralisée (pour des géométries plus complexes) ces grandeurs doivent être calculées et non "encodées en dur" //////////////////////////////////////////
    
    const scalar dt_ = 5e-6;   // pas de temps utilisé pour la simulation, doit être plus petit que le temps de relaxation d'une particule d'une taille donnée
    const scalar UTau = 0.32;		// vitesse de cisaillement pariétale (ici constante en raison de la géométrie très simple) 
    const scalar nu_ = 0.000014607;    // viscosité dynamique
    const scalar dp = 5e-06;		//diamètre des particules considérées

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////   
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    



    Random& rnd = this->owner().rndGen();

    const vector& gradSigma = this->gradSigmaPtr_->primitiveField()[celli];
    const vector& gradSigmaCL = this->gradSigmaCLPtr_->primitiveField()[celli];
    const vector& gradSigmaCL2 = this->gradSigmaCL4Ptr_->primitiveField()[celli];


    scalar y_(wallDist::New(this->owner().mesh()).y()[celli]);

    const scalar k = this->kPtr_->primitiveField()[celli];
    const scalar& epsilon = this->epsilonPtr_->primitiveField()[celli]+rootVSmall;
    
    const symmTensor R = this->sigmaPtr_->primitiveField()[celli];

    const scalar sigmaX = sqrt(R.component(0));
    const scalar sigmaY = sqrt(R.component(3));
    const scalar sigmaZ = sqrt(R.component(5));
    
    const vector RMS_ (sigmaX, sigmaY , sigmaZ);
    
    RMS = RMS_; 
    
    const scalar UV_ = R.component(1);
    
    const scalar gradSigmaX = gradSigma.component(0);
    const scalar gradSigmaY = gradSigma.component(1);
    const scalar gradSigmaZ = gradSigma.component(2);
    
    const scalar gradSigmaCLY = gradSigmaCL.component(1);
    
    const scalar gradSigmaCL2Y = gradSigmaCL2.component(1);   




/////////////////////Calcul de la distance adimensionnée à la paroi de la cellule considérée /////////////////////////////////////////    
            
    scalar yStar = UTau*y_/nu_;  
    
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////Calcul du temps caractéristique lagrangien en fonction de la distance à la paroi //////////////////////////
    
    scalar tL;
    
    if (yStar > 100)
    {
    	tL = k/epsilon;
    	
    }
    else
    {
    	if (yStar > 5)
    	{
    		tL = (5.13+0.27*yStar-0.000474*pow(yStar, 2))*(nu_/pow(UTau, 2)); 
    	}
    	else
    	{
    		tL = 3*(nu_/pow(UTau, 2));  
    	}
    }
    
    
    
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    
    const scalar Kn = (2.0*6.6e-08)/dp; // calcul du nombre de Knudsen associé à la granulométrie des particules, 6.6e-08 (m) représente le libre parcours moyen lambda dans l'air
    const scalar Cc = 1 + Kn*(1.257 + 0.4*exp(-1.1/Kn));     // Facteur de Cunningham pour corriger la force de traînée pour les particules nanométriques    
    const scalar Tp = 2000*pow(dp , 2)*Cc/(18*nu_);     //  Calcul du temps de relaxation de la particule  
    const scalar Stk=Tp*(pow(UTau, 2)/nu_);       //  Calcul du nombre de Stokes pour le terme de correction

    
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////Calcul de 3 nombres aléatoires indépendants pour le terme stochastique ////////////////////////////////////    
    
    
    scalar lambda1 = rnd.scalarNormal();
    scalar lambda2 = rnd.scalarNormal();
    scalar lambda3 = rnd.scalarNormal();


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    const scalar Ki = exp(-dt_/tL);  // tL est supposé être égal dans toutes les directions
    const scalar Ruv = UV_/(sigmaX*sigmaY);
    const scalar b = (Ruv*(1-pow(Ki,2)))/(1-exp(-2*dt_/tL)); // a modifier si on distingue les tL selon la direction
    

////////////////////////////// Calcul du premier terme /////////////////////////////////////////////////////////////////////////////    
    
    
    vector terme1;
    
    if (RMSn.component(0)==0 || RMSn.component(1)==0 || RMSn.component(2)==0)
    {
    	terme1 =  Ki*UTurbn;
    }
    else
    {
    	vector RMSRatio (RMS.component(0)/(RMSn.component(0)+10e-12), RMS.component(1)/(RMSn.component(1)+10e-12), RMS.component(2)/(RMSn.component(2)+10e-12));
    	
    	terme1 = Ki * cmptMultiply(RMSRatio , UTurbn);
    }
    
    
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
////////////////////////////// Calcul du second terme (terme stochastique) /////////////////////////////////////////////////////////    
    
    
    const scalar terme2_x = sigmaX * sqrt(1-pow(Ki,2)) * (sqrt(1-pow(b,2))*lambda1 + b*lambda2);
    const scalar terme2_y = sigmaY * sqrt(1-pow(Ki,2)) * lambda2;
    const scalar terme2_z = sigmaZ * sqrt(1-pow(Ki,2)) * lambda3;
    
    const vector terme2 (terme2_x, terme2_y, terme2_z);
    

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
////////////////////////////// Calcul du terme de correction dans la couche limite  (y+ < 100)   ////////////////////////////////////    
    
       
    const scalar correction_x = RMS.component(0) * (tL/(1+Stk)) * gradSigmaCLY * (1-Ki);
    const scalar correction_y = (tL/(1+Stk))* sigmaY * gradSigmaCL2Y * (1-Ki);
    
    const vector terme3CL (correction_x, correction_y, 0);

    
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////     
////////////////////////////// Calcul du terme de correction hors de la couche limite (y+ > 100) ////////////////////////////////////        
    
    vector terme3;
    
    const scalar terme3_x = (tL/(3*(1+Stk)))* gradSigmaX * (1-Ki);
    const scalar terme3_y = (tL/(3*(1+Stk)))* gradSigmaY * (1-Ki);
    const scalar terme3_z = (tL/(3*(1+Stk)))* gradSigmaZ * (1-Ki);
    
    const vector A (terme3_x, terme3_y, terme3_z);
    
    if (RMSn.component(0)==0 || RMSn.component(1)==0 || RMSn.component(2)==0)
    {
    	terme3 =  A ;
    }
    else
    {
    	vector RMSRatio (RMS.component(0)/(RMSn.component(0)+10e-12), RMS.component(1)/(RMSn.component(1)+10e-12), RMS.component(2)/(RMSn.component(2)+10e-12));
    	
    	terme3 = cmptMultiply(RMSRatio , A);
    }

 
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////    
///////////////// Calcul de la vitesse turbulente selon la distance à la paroi de la cellule ////////////////////////////////////////   
    
    
    if (yStar < 500)
    {
 	UTurb = terme1 + terme2 + terme3CL;
    }
     
    else
    {
     	UTurb = terme1 + terme2 + terme3;
    }


/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


    RMSn   = RMS;      //sauvegarde des valeurs RMS des fluctuations de la vitesse pour leur utilisation à l'itération suivante
    UTurbn = UTurb;	//sauvegarde des valeurs de la vitesse turbulente pour l'itération suivante

    return Uc + UTurb;
}




    
// ************************************************************************* //
