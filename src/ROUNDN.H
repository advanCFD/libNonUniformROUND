/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    Foam::ROUNDNLimiter

Description
    The ROUND formulation on Non-unifrom girds (ROUNDN). 
    Class with limiter function which returns the limiter for the
    ROUNDN differencing scheme based on r obtained from the LimiterFunc
    class.

    Used in conjunction with the template class LimitedScheme.

SourceFiles
    ROUNDN.C

\*---------------------------------------------------------------------------*/

#ifndef ROUNDN_H
#define ROUNDN_H

#include "vector.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ROUNDNLimiter Declaration
\*---------------------------------------------------------------------------*/

template<class LimiterFunc>
class ROUNDNLimiter
:
    public LimiterFunc
{

public:

    ROUNDNLimiter(Istream&)
    {}

    scalar limiter
    (
        const scalar cdWeight,
        const scalar faceFlux,
        const typename LimiterFunc::phiType& phiP,
        const typename LimiterFunc::phiType& phiN,
        const typename LimiterFunc::gradPhiType& gradcP,
        const typename LimiterFunc::gradPhiType& gradcN,
        const vector& d
    ) const
    {
        scalar r = LimiterFunc::r
        (
            faceFlux, phiP, phiN, gradcP, gradcN, d
        );

        scalar s;
        
        if (faceFlux > 0)
    	{         
         s=(1.0-cdWeight)/cdWeight;
    	}
    	else
    	{	 
	 s=cdWeight/(1.0-cdWeight);
    	}

        scalar p2=2.0*(1.0+s+s*r)/(1.0+s)/(1.0+2.0*s);

        scalar rr = 0.0;
        
        if(r>=0.0){
        
        scalar a=(4.0+s)*r/(1.0+s);
        scalar c=(10.0+8.0*s)/(5.0+5.0*s);
        
        rr=min(min(a,p2),c);
        
        }else if(r>-1.0 && r<0.0){
        
        rr=min(-1.0*r,p2);
        
        }else{
        
        rr=min(-0.4,p2);
        
        }
 
        if(faceFlux > 0)
        {
        return rr/(1.0-cdWeight)/2.0;
        }
        else{
        return rr/cdWeight/2.0;
        }
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
