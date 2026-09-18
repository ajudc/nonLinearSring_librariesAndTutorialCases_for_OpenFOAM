/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

InNamespace
    Foam

Description
    Interpolates y values from one curve to another with a different x
    distribution.

    Uses Catmull-Rom spline interpolation between points.

SourceFiles
    interpolateSplineXY.C

\*---------------------------------------------------------------------------*/

#ifndef interExtrapolateSplineXY_H
#define interExtrapolateSplineXY_H

#include "Field.H"
#include "scalarField.H"
#include "Enum.H"

namespace Foam
{

//- Extrapolation behaviour outside bounds
enum splineExtrapolationType
{
    SET_linear,
    SET_cubic
};

//- Name ↔ enum mapping
extern const Enum<splineExtrapolationType>
    splineExtrapolationTypeNames;


//- Interpolate a field of values
template<class Type>
Field<Type> interExtrapolateSplineXY
(
    const scalarField& xNew,
    const scalarField& xOld,
    const Field<Type>& yOld,
    const splineExtrapolationType extrapType = SET_cubic
);


//- Interpolate a single value
template<class Type>
Type interExtrapolateSplineXY
(
    const scalar x,
    const scalarField& xOld,
    const Field<Type>& yOld,
    const splineExtrapolationType extrapType = SET_cubic
);

} // End namespace Foam

#endif
