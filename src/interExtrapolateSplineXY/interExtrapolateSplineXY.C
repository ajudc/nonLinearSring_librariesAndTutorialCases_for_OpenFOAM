/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011 OpenFOAM Foundation
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

#include "interExtrapolateSplineXY.H"
#include "primitiveFields.H"
#include "error.H"

namespace Foam
{

// * * * * * * * * * * * * * Enum definitions * * * * * * * * * * * * * //

const Enum<splineExtrapolationType>
splineExtrapolationTypeNames
({
    { SET_linear, "linearExtrapolation" },
    { SET_cubic,  "cubicExtrapolation" }
});


// * * * * * * * * * * * * * Field version * * * * * * * * * * * * * //

template<class Type>
Field<Type> interExtrapolateSplineXY
(
    const scalarField& xNew,
    const scalarField& xOld,
    const Field<Type>& yOld,
    const splineExtrapolationType extrapType
)
{
    Field<Type> yNew(xNew.size());

    forAll(xNew, i)
    {
        yNew[i] =
            interExtrapolateSplineXY
            (
                xNew[i],
                xOld,
                yOld,
                extrapType
            );
    }

    return yNew;
}


// * * * * * * * * * * * * * Scalar version * * * * * * * * * * * * * //

template<class Type>
Type interExtrapolateSplineXY
(
    const scalar x,
    const scalarField& xOld,
    const Field<Type>& yOld,
    const splineExtrapolationType extrapType
)
{
    const label n = xOld.size();

    if (n == 0)
    {
        FatalErrorInFunction
            << "Cannot interpolate from an empty table"
            << exit(FatalError);
    }

    if (n == 1)
    {
        return yOld[0];
    }

    // ---------- LINEAR EXTRAPOLATION ----------
    if (extrapType == SET_linear)
    {
        if (x <= xOld[0])
        {
            return
                yOld[0]
              + (x - xOld[0])
               *(yOld[1] - yOld[0])
               /(xOld[1] - xOld[0]);
        }

        if (x >= xOld[n-1])
        {
            return
                yOld[n-1]
              + (x - xOld[n-1])
               *(yOld[n-1] - yOld[n-2])
               /(xOld[n-1] - xOld[n-2]);
        }
    }

    // ---------- CUBIC INTERPOLATION / EXTRAPOLATION ----------

    label lo = -1;
    label hi = -1;

    if (x <= xOld[0])
    {
        lo = 0;
        hi = 1;
    }
    else if (x >= xOld[n-1])
    {
        lo = n - 2;
        hi = n - 1;
    }
    else
    {
        for (label i = 1; i < n; ++i)
        {
            if (xOld[i] >= x)
            {
                lo = i - 1;
                hi = i;
                break;
            }
        }
    }

    const Type& y1 = yOld[lo];
    const Type& y2 = yOld[hi];

    Type y0 = (lo == 0 ? 2*y1 - y2 : yOld[lo - 1]);
    Type y3 = (hi + 1 >= n ? 2*y2 - y1 : yOld[hi + 1]);

    const scalar mu =
        (x - xOld[lo])/(xOld[hi] - xOld[lo]);

    return
        0.5
       *(
            2*y1
          + mu
           *(
               -y0 + y2
              + mu
               *(
                   (2*y0 - 5*y1 + 4*y2 - y3)
                 + mu*(-y0 + 3*y1 - 3*y2 + y3)
               )
            )
        );
}

// * * * * * * * * * * * * Explicit template instantiation * * * * * * * * * //
template scalar interExtrapolateSplineXY<scalar>
(
    const scalar,
    const scalarField&,
    const Field<scalar>&,
    const splineExtrapolationType
);

template Field<scalar> interExtrapolateSplineXY<scalar>
(
    const scalarField&,
    const scalarField&,
    const Field<scalar>&,
    const splineExtrapolationType
);


} // End namespace Foam
