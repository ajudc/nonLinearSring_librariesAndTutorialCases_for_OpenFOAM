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
    Foam::sixDoFRigidBodyMotionRestraints::nonLinearSpring

Description
    sixDoFRigidBodyMotionRestraints model. Nonlinear spring with damping.

SourceFiles
    nonLinearSpring.C

\*---------------------------------------------------------------------------*/

#ifndef nonLinearSpring_H
#define nonLinearSpring_H

#include "sixDoFRigidBodyMotionRestraint.H"
#include "interExtrapolationTable.H"
#include "interExtrapolateSplineXY.H"
#include "autoPtr.H"
#include "vector.H"
#include "scalarField.H"

namespace Foam
{
namespace sixDoFRigidBodyMotionRestraints
{

class nonLinearSpring
:
    public sixDoFRigidBodyMotionRestraint
{
    // Private data

        vector anchor_;
        vector refAttachmentPt_;
        scalar restLength_;

        // Spring model
        word modelType_;
        word stiffnessInterpolationType_;
        scalar A_;
        scalar B_;

        // Damping model
        word dampingModelType_;
        word dampingInterpolationType_;
        word listType_;
        scalar c1_;
        scalar c2_;

        // Linear tables
        autoPtr<interExtrapolationTable<scalar>> stiffnessTablePtr_;
        autoPtr<interExtrapolationTable<scalar>> dampingTablePtr_;

        // Cubic spline storage
        autoPtr<scalarField> springX_;
        autoPtr<scalarField> springY_;
        autoPtr<scalarField> dampingX_;
        autoPtr<scalarField> dampingY_;

        splineExtrapolationType springExtrapType_;
        splineExtrapolationType dampingExtrapType_;


public:

    TypeName("nonLinearSpring");

    nonLinearSpring
    (
        const word& name,
        const dictionary& dict
    );

    virtual ~nonLinearSpring();

    virtual autoPtr<sixDoFRigidBodyMotionRestraint> clone() const;

    virtual void restrain
    (
        const sixDoFRigidBodyMotion& motion,
        vector& restraintPosition,
        vector& restraintForce,
        vector& restraintMoment
    ) const;

    virtual bool read(const dictionary& dict);

    virtual void write(Ostream& os) const;
};

} // End namespace sixDoFRigidBodyMotionRestraints
} // End namespace Foam

#endif


