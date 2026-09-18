/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Custom interpolation/extrapolation table based on interpolationTable
-------------------------------------------------------------------------------
License
    This file is not part of the official OpenFOAM distribution.

    It is derived from interpolationTable and extended to support an
    additional bounding mode (EXTRAPOLATE).

\*---------------------------------------------------------------------------*/

#ifndef interExtrapolationTable_H
#define interExtrapolationTable_H

#include "List.H"
#include "Tuple2.H"
#include "tableReader.H"
#include "autoPtr.H"
#include "Field.H"
#include "fileName.H"
#include "UList.H"
#include "Enum.H"     

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Bounding modes for interExtrapolationTable
\*---------------------------------------------------------------------------*/

//- Handling for out-of-bound values for interExtrapolationTable
enum interExtrapolationBounding
{
    IEB_ERROR,
    IEB_WARN,
    IEB_CLAMP,
    IEB_REPEAT,
    IEB_EXTRAPOLATE
};

//- Name ↔ enum mapping for interExtrapolationBounding (for dictionary I/O)
extern const Enum<interExtrapolationBounding>
    interExtrapolationBoundingNames;


/*---------------------------------------------------------------------------*\
                Class interExtrapolationTable Declaration
\*---------------------------------------------------------------------------*/

template<class Type>
class interExtrapolationTable
:
    public List<Tuple2<scalar, Type>>
{
    // Private Data

        //- Handling for out-of-bound values
        interExtrapolationBounding bounding_;

        //- File name
        fileName fileName_;

        //- Table reader
        autoPtr<tableReader<Type>> reader_;


    // Private Member Functions

        //- Read the table of data from file
        void readTable();


public:

    // Public Data Types

        //- The element data type
        typedef Tuple2<scalar, Type> value_type;

        //- The mapped data type
        typedef Type mapped_type;


    // Constructors

        //- Default construct
        interExtrapolationTable();

        //- Construct from components
        interExtrapolationTable
        (
            const List<Tuple2<scalar, Type>>& values,
            const interExtrapolationBounding bounding,
            const fileName& fName
        );

        //- Construct given the name of the file containing the table of data
        explicit interExtrapolationTable(const fileName& fName);

        //- Construct by reading file name and outOfBounds from dictionary
        //- and read the table from that file.
        //  This is a specialised constructor used by patchFields
        explicit interExtrapolationTable(const dictionary& dict);

        //- Copy construct
        interExtrapolationTable(const interExtrapolationTable& tbl);


    // Member Functions

        //- Return an interpolated value in List
        static Type interpolateValue
        (
            const List<Tuple2<scalar, Type>>& list,
            scalar lookupValue,
            interExtrapolationBounding bounding = IEB_CLAMP
        );

        //- Check that list is monotonically increasing
        //  Exit with a FatalError if there is a problem
        void check() const;

        //- Write
        void write(Ostream& os) const;

        //- Return the rate of change at the interpolation location
        //- for the given lookup value
        Type rateOfChange(scalar lookupValue) const;

        //- Return an interpolated value
        Type interpolateValue(scalar lookupValue) const;

        //- Return multiple interpolated values
        tmp<Field<Type>> interpolateValues
        (
            const UList<scalar>& vals
        ) const;


    // Member Operators

        //- Copy assignment
        void operator=(const interExtrapolationTable<Type>& rhs);

        //- Return an element of constant Tuple2<scalar, Type>
        const Tuple2<scalar, Type>& operator[](label idx) const;

        //- Return an interpolated value
        Type operator()(scalar lookupValue) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

#endif

// ************************************************************************* //
