/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2018-2020 OpenCFD Ltd.
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

#include "nonLinearAxialAngularSpring.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotion.H"
#include "transform.H"

namespace Foam
{
namespace sixDoFRigidBodyMotionRestraints
{

defineTypeNameAndDebug(nonLinearAxialAngularSpring, 0);

addToRunTimeSelectionTable
(
    sixDoFRigidBodyMotionRestraint,
    nonLinearAxialAngularSpring,
    dictionary
);

}
}


// * * * * * * * * * * * * * Constructors * * * * * * * * * * * * * //

Foam::sixDoFRigidBodyMotionRestraints::nonLinearAxialAngularSpring::
nonLinearAxialAngularSpring
(
    const word& name,
    const dictionary& dict
)
:
    sixDoFRigidBodyMotionRestraint(name, dict),
    refQ_(I),
    axis_(Zero),
    modelType_("stiffnessList"),
    A_(0.0),
    B_(0.0),
    dampingModelType_("linear"),
    listType_("angularVelocity"),
    c1_(0.0),
    c2_(0.0),
    stiffnessInterpolationType_("linear"),
    dampingInterpolationType_("linear"),
    stiffnessTable_(nullptr),
    dampingTable_(nullptr)
{
    read(dict);
}


Foam::sixDoFRigidBodyMotionRestraints::nonLinearAxialAngularSpring::
~nonLinearAxialAngularSpring() = default;


Foam::autoPtr<Foam::sixDoFRigidBodyMotionRestraint>
Foam::sixDoFRigidBodyMotionRestraints::nonLinearAxialAngularSpring::clone() const
{
    return autoPtr<sixDoFRigidBodyMotionRestraint>
    (
        new nonLinearAxialAngularSpring(*this)
    );
}


// * * * * * * * * * * * * * Restrain * * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotionRestraints::nonLinearAxialAngularSpring::restrain
(
    const sixDoFRigidBodyMotion& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    vector refDir = rotationTensor(vector(1,0,0), axis_) & vector(0,1,0);

    vector oldDir = refQ_ & refDir;
    vector newDir = motion.orientation() & refDir;

    oldDir -= (axis_ & oldDir)*axis_;
    newDir -= (axis_ & newDir)*axis_;

    oldDir /= mag(oldDir) + VSMALL;
    newDir /= mag(newDir) + VSMALL;

    scalar theta = acos(min(oldDir & newDir, 1.0));

    vector a = (oldDir ^ newDir);
    a = (a & axis_)*axis_;
    a /= mag(a) + VSMALL;

    // -------- SPRING --------
    scalar springMag = 0.0;

    if (modelType_ == "linear")
    {
        springMag = A_*theta;
    }
    else if (modelType_ == "Duffing")
    {
        springMag = A_*theta + B_*pow3(theta);
    }
    else if (modelType_ == "stiffnessList")
    {
        if (stiffnessInterpolationType_ == "linear")
        {
            springMag = stiffnessTable_()(theta);
        }
        else
        {
            springMag =
                interExtrapolateSplineXY
                (
                    theta,
                    stiffnessX_,
                    stiffnessY_,
                    SET_cubic
                );
        }
    }

    // -------- DAMPING --------
    scalar omega = (motion.omega() & a);
    scalar input =
        (listType_ == "angularVelocity" ? omega : theta);

    scalar dampingMag = 0.0;

    if (dampingModelType_ == "linear")
    {
        dampingMag = c1_*omega;
    }
    else if (dampingModelType_ == "Rayleigh")
    {
        dampingMag = c1_*omega + c2_*pow3(omega);
    }
    else if (dampingModelType_ == "Mann")
    {
        dampingMag = c1_*omega + c2_*sqr(theta)*omega;
    }
    else if (dampingModelType_ == "dampingList")
    {
        if (dampingInterpolationType_ == "linear")
        {
            dampingMag = dampingTable_()(input);
        }
        else
        {
            dampingMag =
                interExtrapolateSplineXY
                (
                    input,
                    dampingX_,
                    dampingY_,
                    SET_cubic
                );
        }
    }

    restraintMoment = -(springMag + dampingMag)*a;
    restraintForce  = Zero;
    restraintPosition = motion.centreOfRotation();
}


// * * * * * * * * * * * * * Read dictionary * * * * * * * * * * * * * //

bool Foam::sixDoFRigidBodyMotionRestraints::nonLinearAxialAngularSpring::read
(
    const dictionary& dict
)
{
    sixDoFRigidBodyMotionRestraint::read(dict);

    refQ_ = sDoFRBMRCoeffs_.getOrDefault<tensor>("referenceOrientation", I);
    sDoFRBMRCoeffs_.readEntry("axis", axis_);
    axis_ /= mag(axis_) + VSMALL;

    modelType_ =
        sDoFRBMRCoeffs_.lookupOrDefault<word>("modelType", "stiffnessList");

    if (modelType_ == "linear")
    {
        sDoFRBMRCoeffs_.readEntry("A", A_);
    }
    else if (modelType_ == "Duffing")
    {
        sDoFRBMRCoeffs_.readEntry("A", A_);
        sDoFRBMRCoeffs_.readEntry("B", B_);
    }
    else if (modelType_ == "stiffnessList")
    {
        stiffnessInterpolationType_ =
            sDoFRBMRCoeffs_.lookupOrDefault<word>
            (
                "interpolationType",
                "linear"
            );

        List<Tuple2<scalar, scalar>> values;
        sDoFRBMRCoeffs_.subDict("stiffnessList").lookup("values") >> values;

        if (stiffnessInterpolationType_ == "linear")
        {
            stiffnessTable_.reset
            (
                new interExtrapolationTable<scalar>
                (
                    values,
                    IEB_EXTRAPOLATE,
                    "nonlinearAxialAngularSpring::stiffness"
                )
            );
        }
        else
        {
            stiffnessX_.setSize(values.size());
            stiffnessY_.setSize(values.size());

            forAll(values, i)
            {
                stiffnessX_[i] = values[i].first();
                stiffnessY_[i] = values[i].second();
            }
        }
    }

    dampingModelType_ =
        sDoFRBMRCoeffs_.lookupOrDefault<word>
        (
            "dampingModelType",
            "linear"
        );

    if (dampingModelType_ == "linear")
    {
        sDoFRBMRCoeffs_.readEntry("c1", c1_);
    }
    else if (dampingModelType_ == "Rayleigh" || dampingModelType_ == "Mann")
    {
        sDoFRBMRCoeffs_.readEntry("c1", c1_);
        sDoFRBMRCoeffs_.readEntry("c2", c2_);
    }
    else if (dampingModelType_ == "dampingList")
    {
        sDoFRBMRCoeffs_.readEntry("listType", listType_);

        dampingInterpolationType_ =
            sDoFRBMRCoeffs_.lookupOrDefault<word>
            (
                "dampingInterpolationType",
                "linear"
            );

        List<Tuple2<scalar, scalar>> values;
        sDoFRBMRCoeffs_.subDict("dampingList").lookup("values") >> values;

        if (dampingInterpolationType_ == "linear")
        {
            dampingTable_.reset
            (
                new interExtrapolationTable<scalar>
                (
                    values,
                    IEB_EXTRAPOLATE,
                    "nonLinearAxialAngularSpring::damping"
                )
            );
        }
        else
        {
            dampingX_.setSize(values.size());
            dampingY_.setSize(values.size());

            forAll(values, i)
            {
                dampingX_[i] = values[i].first();
                dampingY_[i] = values[i].second();
            }
        }
    }

    return true;
}


// * * * * * * * * * * * * * Write * * * * * * * * * * * * * //

void Foam::sixDoFRigidBodyMotionRestraints::nonLinearAxialAngularSpring::write
(
    Ostream& os
) const
{
    os.writeEntry("referenceOrientation", refQ_);
    os.writeEntry("axis", axis_);
    os.writeEntry("modelType", modelType_);
    os.writeEntry("dampingModelType", dampingModelType_);
}




