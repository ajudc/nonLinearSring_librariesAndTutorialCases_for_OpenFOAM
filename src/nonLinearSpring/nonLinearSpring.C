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

\*---------------------------------------------------------------------------*/

#include "nonLinearSpring.H"
#include "addToRunTimeSelectionTable.H"
#include "sixDoFRigidBodyMotion.H"
#include "dictionary.H"

namespace Foam
{
namespace sixDoFRigidBodyMotionRestraints
{

    defineTypeNameAndDebug(nonLinearSpring, 0);
    
    addToRunTimeSelectionTable
    (
        sixDoFRigidBodyMotionRestraint,
        nonLinearSpring,
        dictionary
    );


// * * * * * * * * * * * * Constructors * * * * * * * * * * * * //

nonLinearSpring::nonLinearSpring
(
    const word& name,
    const dictionary& dict
)
:
    sixDoFRigidBodyMotionRestraint(name, dict),
    anchor_(Zero),
    refAttachmentPt_(Zero),
    restLength_(0.0),
    modelType_("stiffnessList"),
    stiffnessInterpolationType_("linear"),
    A_(0.0),
    B_(0.0),
    dampingModelType_("linear"),
    dampingInterpolationType_("linear"),
    listType_("velocity"),
    c1_(0.0),
    c2_(0.0),
    stiffnessTablePtr_(nullptr),
    dampingTablePtr_(nullptr),
    springExtrapType_(SET_cubic),
    dampingExtrapType_(SET_cubic)
{
    read(dict);
}

nonLinearSpring::~nonLinearSpring() = default;


autoPtr<sixDoFRigidBodyMotionRestraint>
nonLinearSpring::clone() const
{
    return autoPtr<sixDoFRigidBodyMotionRestraint>
    (
        new nonLinearSpring(*this)
    );
}


// * * * * * * * * * * * * Restraint * * * * * * * * * * * * //

void nonLinearSpring::restrain
(
    const sixDoFRigidBodyMotion& motion,
    vector& restraintPosition,
    vector& restraintForce,
    vector& restraintMoment
) const
{
    restraintPosition = motion.transform(refAttachmentPt_);

    vector r = restraintPosition - anchor_;
    scalar magR = mag(r);
    r /= (magR + VSMALL);

    vector v = motion.velocity(restraintPosition);
    scalar displacement = magR - restLength_;

    scalar springForce = 0.0;

    // ---------- SPRING ----------
    if (modelType_ == "stiffnessList")
    {
        if (stiffnessInterpolationType_ == "linear")
        {
            springForce = stiffnessTablePtr_()(displacement);
        }
        else
        {
            springForce =
                interExtrapolateSplineXY
                (
                    displacement,
                    *springX_,
                    *springY_,
                    springExtrapType_
                );
        }
    }
    else if (modelType_ == "linear")
    {
        springForce = A_*displacement;
    }
    else if (modelType_ == "Duffing")
    {
        springForce = A_*displacement + B_*pow3(displacement);
    }
    else
    {
        FatalErrorInFunction
            << "Unknown modelType: " << modelType_
            << abort(FatalError);
    }

    // ---------- DAMPING ----------
    scalar vel = (r & v);
    scalar dampingForce = 0.0;

    if (dampingModelType_ == "linear")
    {
        dampingForce = c1_*vel;
    }
    else if (dampingModelType_ == "Rayleigh")
    {
        dampingForce = c1_*vel + c2_*pow3(vel);
    }
    else if (dampingModelType_ == "Mann")
    {
        dampingForce = c1_*vel + c2_*sqr(displacement)*vel;
    }
    else if (dampingModelType_ == "dampingList")
    {
        scalar input =
            (listType_ == "velocity" ? vel : displacement);

        if (dampingInterpolationType_ == "linear")
        {
            dampingForce = (*dampingTablePtr_)(input);
        }
        else
        {
            dampingForce =
                interExtrapolateSplineXY
                (
                    input,
                    *dampingX_,
                    *dampingY_,
                    dampingExtrapType_
                );
        }
    }
    else
    {
        FatalErrorInFunction
            << "Unknown dampingModelType: " << dampingModelType_
            << abort(FatalError);
    }

    restraintForce = -(springForce + dampingForce)*r;
    restraintMoment = Zero;
}


// * * * * * * * * * * * * Read dictionary * * * * * * * * * * * //

bool nonLinearSpring::read(const dictionary& dict)
{
    sixDoFRigidBodyMotionRestraint::read(dict);

    sDoFRBMRCoeffs_.readEntry("anchor", anchor_);
    sDoFRBMRCoeffs_.readEntry("refAttachmentPt", refAttachmentPt_);
    sDoFRBMRCoeffs_.readEntry("restLength", restLength_);

    modelType_ =
        sDoFRBMRCoeffs_.lookupOrDefault<word>
        (
            "modelType",
            "stiffnessList"
        );

    // ---------- SPRING ----------
    if (modelType_ == "stiffnessList")
    {
        const dictionary& kDict =
            sDoFRBMRCoeffs_.subDict("stiffnessList");

        stiffnessInterpolationType_ =
            kDict.lookupOrDefault<word>
            (
                "interpolationType",
                "linear"
            );

        List<Tuple2<scalar, scalar>> values;
        kDict.lookup("values") >> values;

        if (stiffnessInterpolationType_ == "linear")
        {
            interExtrapolationBounding bound =
                interExtrapolationBoundingNames.getOrDefault
                (
                    "outOfBounds",
                    kDict,
                    IEB_EXTRAPOLATE,
                    true
                );

            stiffnessTablePtr_.reset
            (
                new interExtrapolationTable<scalar>
                (
                    values,
                    bound,
                    "nonLinearSpring::stiffness"
                )
            );
            springX_.clear();
            springY_.clear();
        }
        else
        {
            springX_.reset(new scalarField(values.size()));
            springY_.reset(new scalarField(values.size()));

            forAll(values, i)
            {
                (*springX_)[i] = values[i].first();
                (*springY_)[i] = values[i].second();
            }

            word extrapName =
                kDict.lookupOrDefault<word>
                (
                    "outOfBounds",
                    "cubicExtrapolation"
                );

            springExtrapType_ =
                splineExtrapolationTypeNames[extrapName];

            stiffnessTablePtr_.clear();
        }
    }
    else if (modelType_ == "linear")
    {
        sDoFRBMRCoeffs_.readEntry("A", A_);
    }
    else if (modelType_ == "Duffing")
    {
        sDoFRBMRCoeffs_.readEntry("A", A_);
        sDoFRBMRCoeffs_.readEntry("B", B_);
    }

    // ---------- DAMPING ----------
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
    else if (dampingModelType_ == "Rayleigh"
          || dampingModelType_ == "Mann")
    {
        sDoFRBMRCoeffs_.readEntry("c1", c1_);
        sDoFRBMRCoeffs_.readEntry("c2", c2_);
    }
    else if (dampingModelType_ == "dampingList")
    {
        sDoFRBMRCoeffs_.readEntry("listType", listType_);

        const dictionary& dDict =
            sDoFRBMRCoeffs_.subDict("dampingList");

        dampingInterpolationType_ =
            dDict.lookupOrDefault<word>
            (
                "interpolationType",
                "linear"
            );

        List<Tuple2<scalar, scalar>> values;
        dDict.lookup("values") >> values;

        if (dampingInterpolationType_ == "linear")
        {
            interExtrapolationBounding bound =
                interExtrapolationBoundingNames.getOrDefault
                (
                    "outOfBounds",
                    dDict,
                    IEB_EXTRAPOLATE,
                    true
                );

            dampingTablePtr_.reset
            (
                new interExtrapolationTable<scalar>
                (
                    values,
                    bound,
                    "nonLinearSpring::damping"
                )
            );
            dampingX_.clear();
            dampingY_.clear();
        }
        else
        {
            dampingX_.reset(new scalarField(values.size()));
            dampingY_.reset(new scalarField(values.size()));

            forAll(values, i)
            {
                (*dampingX_)[i] = values[i].first();
                (*dampingY_)[i] = values[i].second();
            }

            word extrapName =
                dDict.lookupOrDefault<word>
                (
                    "outOfBounds",
                    "cubicExtrapolation"
                );

            dampingExtrapType_ =
                splineExtrapolationTypeNames[extrapName];

            dampingTablePtr_.clear();
        }
    }

    Info<< "nonLinearSpring loaded: modelType=" << modelType_
        << ", dampingModelType=" << dampingModelType_ << endl;

    return true;
}


void nonLinearSpring::write(Ostream& os) const
{
    os.writeEntry("anchor", anchor_);
    os.writeEntry("refAttachmentPt", refAttachmentPt_);
    os.writeEntry("restLength", restLength_);
    os.writeEntry("modelType", modelType_);
    os.writeEntry("dampingModelType", dampingModelType_);
}

} // namespace sixDoFRigidBodyMotionRestraints
} // namespace Foam

