/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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
    Foam::sixDoFRigidBodyMotionRestraints::linearAxialAngularSpring

Description
    sixDoFRigidBodyMotionRestraints model.  Linear axial angular spring.

SourceFiles
    linearAxialAngularSpring.C

\*---------------------------------------------------------------------------*/

#ifndef nonLinearAxialAngularSpring_H
#define nonLinearAxialAngularSpring_H

#include "sixDoFRigidBodyMotionRestraint.H"
#include "interExtrapolationTable.H"
#include "interExtrapolateSplineXY.H"
#include "vector.H"
#include "tensor.H"

namespace Foam
{
namespace sixDoFRigidBodyMotionRestraints
{

class nonLinearAxialAngularSpring
:
    public sixDoFRigidBodyMotionRestraint
{
    // Reference orientation
    tensor refQ_;

    // Rotation axis
    vector axis_;

    // Spring model
    word modelType_;
    scalar A_;
    scalar B_;

    // Damping model
    word dampingModelType_;
    word listType_;
    scalar c1_;
    scalar c2_;

    // Interpolation control
    word stiffnessInterpolationType_;
    word dampingInterpolationType_;

    // Linear interpolation tables
    autoPtr<interExtrapolationTable<scalar>> stiffnessTable_;
    autoPtr<interExtrapolationTable<scalar>> dampingTable_;

    // Cubic spline data
    scalarField stiffnessX_;
    scalarField stiffnessY_;
    scalarField dampingX_;
    scalarField dampingY_;

public:

    TypeName("nonLinearAxialAngularSpring");

    // Constructors
    nonLinearAxialAngularSpring
    (
        const word& name,
        const dictionary& dict
    );

    // Destructor
    virtual ~nonLinearAxialAngularSpring();

    // Clone
    virtual autoPtr<sixDoFRigidBodyMotionRestraint> clone() const;

    // Main restraint calculation
    virtual void restrain
    (
        const sixDoFRigidBodyMotion& motion,
        vector& restraintPosition,
        vector& restraintForce,
        vector& restraintMoment
    ) const;

    // Dictionary I/O
    virtual bool read(const dictionary& dict);
    virtual void write(Ostream& os) const;
};

} // End namespace sixDoFRigidBodyMotionRestraints
} // End namespace Foam

#endif

