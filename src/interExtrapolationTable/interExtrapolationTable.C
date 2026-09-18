/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Implementation of Foam::interExtrapolationTable
    (derived from interpolationTable, extended with EXTRAPOLATE mode)
\*---------------------------------------------------------------------------*/

#include "interExtrapolationTable.H"
#include "openFoamTableReader.H"

namespace Foam
{

// Define names for interExtrapolationBounding (for I/O)
const Enum<interExtrapolationBounding>
    interExtrapolationBoundingNames
(
    {
        { IEB_ERROR,       "ERROR" },
        { IEB_WARN,        "WARN" },
        { IEB_CLAMP,       "CLAMP" },
        { IEB_REPEAT,      "REPEAT" },
        { IEB_EXTRAPOLATE, "EXTRAPOLATE" }
    }
);


// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

template<class Type>
void interExtrapolationTable<Type>::readTable()
{
    // preserve the original (unexpanded) fileName to avoid absolute paths
    // appearing subsequently in the write() method
    fileName fName(fileName_);

    fName.expand();

    // Read data from file
    reader_()(fName, *this);

    if (this->empty())
    {
        FatalErrorInFunction
            << "table read from " << fName << " is empty" << nl
            << exit(FatalError);
    }

    // Check that the data are okay
    check();
}


// * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * * * * //

template<class Type>
interExtrapolationTable<Type>::interExtrapolationTable()
:
    List<value_type>(),
    bounding_(IEB_WARN),
    fileName_("fileNameIsUndefined"),
    reader_(nullptr)
{}


template<class Type>
interExtrapolationTable<Type>::interExtrapolationTable
(
    const List<Tuple2<scalar, Type>>& values,
    const interExtrapolationBounding bounding,
    const fileName& fName
)
:
    List<value_type>(values),
    bounding_(bounding),
    fileName_(fName),
    reader_(nullptr)
{}


template<class Type>
interExtrapolationTable<Type>::interExtrapolationTable(const fileName& fName)
:
    List<value_type>(),
    bounding_(IEB_WARN),
    fileName_(fName),
    reader_(new openFoamTableReader<Type>())
{
    readTable();
}


template<class Type>
interExtrapolationTable<Type>::interExtrapolationTable(const dictionary& dict)
:
    List<value_type>(),
    bounding_
    (
        interExtrapolationBoundingNames.getOrDefault
        (
            "outOfBounds",
            dict,
            IEB_WARN,
            true  // Failsafe behaviour
        )
    ),
    fileName_(dict.get<fileName>("file")),
    reader_(tableReader<Type>::New(dict))
{
    readTable();
}


template<class Type>
interExtrapolationTable<Type>::interExtrapolationTable
(
    const interExtrapolationTable& tbl
)
:
    List<value_type>(tbl),
    bounding_(tbl.bounding_),
    fileName_(tbl.fileName_),
    reader_(tbl.reader_.clone())
{}


// * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * * * //

template<class Type>
void interExtrapolationTable<Type>::check() const
{
    const List<value_type>& list = *this;

    scalar prevValue(0);

    label i = 0;
    for (const auto& item : list)
    {
        const scalar& currValue = item.first();

        // Avoid duplicate values (divide-by-zero error)
        if (i && currValue <= prevValue)
        {
            FatalErrorInFunction
                << "out-of-order value: "
                << currValue << " at index " << i << nl
                << exit(FatalError);
        }
        prevValue = currValue;
        ++i;
    }
}


template<class Type>
void interExtrapolationTable<Type>::write(Ostream& os) const
{
    os.writeEntry("file", fileName_);
    os.writeEntry("outOfBounds", interExtrapolationBoundingNames[bounding_]);
    if (reader_)
    {
        reader_->write(os);
    }
}


template<class Type>
Type interExtrapolationTable<Type>::rateOfChange(scalar lookupValue) const
{
    const List<value_type>& list = *this;

    const label n = list.size();

    if (n <= 1)
    {
        // Not enough entries for a rate of change
        return Zero;
    }

    const scalar minLimit = list.first().first();
    const scalar maxLimit = list.last().first();

    if (lookupValue < minLimit)
    {
        switch (bounding_)
        {
            case IEB_ERROR:
            {
                FatalErrorInFunction
                    << "value (" << lookupValue << ") less than lower "
                    << "bound (" << minLimit << ")\n"
                    << exit(FatalError);
                break;
            }
            case IEB_WARN:
            {
                WarningInFunction
                    << "value (" << lookupValue << ") less than lower "
                    << "bound (" << minLimit << ")\n"
                    << "    Zero rate of change." << endl;

                // Behaviour as per CLAMP
                return Zero;
            }
            case IEB_CLAMP:
            {
                return Zero;
            }
            case IEB_REPEAT:
            {
                // Adjust lookupValue to >= minLimit
                scalar span = maxLimit - minLimit;
                lookupValue = fmod(lookupValue - minLimit, span) + minLimit;
                break;
            }
            case IEB_EXTRAPOLATE:
            {
                // Use slope of first segment
                const scalar x0 = list[0].first();
                const scalar x1 = list[1].first();
                const Type&  y0 = list[0].second();
                const Type&  y1 = list[1].second();

                return (y1 - y0)/(x1 - x0);
            }
        }
    }
    else if (lookupValue >= maxLimit)
    {
        switch (bounding_)
        {
            case IEB_ERROR:
            {
                FatalErrorInFunction
                    << "value (" << lookupValue << ") greater than upper "
                    << "bound (" << maxLimit << ")\n"
                    << exit(FatalError);
                break;
            }
            case IEB_WARN:
            {
                WarningInFunction
                    << "value (" << lookupValue << ") greater than upper "
                    << "bound (" << maxLimit << ")\n"
                    << "    Zero rate of change." << endl;

                // Behaviour as per CLAMP
                return Zero;
            }
            case IEB_CLAMP:
            {
                return Zero;
            }
            case IEB_REPEAT:
            {
                // Adjust lookupValue <= maxLimit
                scalar span = maxLimit - minLimit;
                lookupValue = fmod(lookupValue - minLimit, span) + minLimit;
                break;
            }
            case IEB_EXTRAPOLATE:
            {
                // Use slope of last segment
                const scalar xNm1 = list[n-2].first();
                const scalar xN   = list[n-1].first();
                const Type&  yNm1 = list[n-2].second();
                const Type&  yN   = list[n-1].second();

                return (yN - yNm1)/(xN - xNm1);
            }
        }
    }

    label lo = 0;
    label hi = 0;

    // Look for the correct range
    for (label i = 0; i < n; ++i)
    {
        if (lookupValue >= list[i].first())
        {
            lo = hi = i;
        }
        else
        {
            hi = i;
            break;
        }
    }

    if (lo == hi)
    {
        return Zero;
    }
    else if (hi == 0)
    {
        // This treatment should only occur under these conditions:
        //  -> the 'REPEAT' treatment
        //  -> (0 <= value <= minLimit)
        //  -> minLimit > 0
        // Use the value at maxLimit as the value for value=0
        lo = n - 1;

        return
        (
            (list[hi].second() - list[lo].second())
          / (list[hi].first() + minLimit - list[lo].first())
        );
    }

    // Normal rate of change
    return
    (
        (list[hi].second() - list[lo].second())
      / (list[hi].first() - list[lo].first())
    );
}


template<class Type>
Type interExtrapolationTable<Type>::interpolateValue
(
    const List<Tuple2<scalar, Type>>& list,
    scalar lookupValue,
    interExtrapolationBounding bounding
)
{
    const label n = list.size();

    if (n <= 1)
    {
        #ifdef FULLDEBUG
        if (!n)
        {
            FatalErrorInFunction
                << "Cannot interpolate from zero-sized table" << nl
                << exit(FatalError);
        }
        #endif

        return list.first().second();
    }

    const scalar minLimit = list.first().first();
    const scalar maxLimit = list.last().first();

    if (lookupValue < minLimit)
    {
        switch (bounding)
        {
            case IEB_ERROR:
            {
                FatalErrorInFunction
                    << "value (" << lookupValue << ") less than lower "
                    << "bound (" << minLimit << ")\n"
                    << exit(FatalError);
                break;
            }
            case IEB_WARN:
            {
                WarningInFunction
                    << "value (" << lookupValue << ") less than lower "
                    << "bound (" << minLimit << ")\n"
                    << "    Continuing with the first entry" << endl;

                // Behaviour as per CLAMP
                return list.first().second();
            }
            case IEB_CLAMP:
            {
                return list.first().second();
            }
            case IEB_REPEAT:
            {
                // adjust lookupValue to >= minLimit
                const scalar span = maxLimit - minLimit;
                lookupValue = fmod(lookupValue - minLimit, span) + minLimit;
                break;
            }
            case IEB_EXTRAPOLATE:
            {
                // Linear extrapolation using first two points
                const scalar x0 = list[0].first();
                const scalar x1 = list[1].first();
                const Type&  y0 = list[0].second();
                const Type&  y1 = list[1].second();

                const scalar dx  = x1 - x0;
                const scalar dxi = 1.0/dx;

                return
                (
                    y0
                  + (y1 - y0)*(lookupValue - x0)*dxi
                );
            }
        }
    }
    else if (lookupValue >= maxLimit)
    {
        switch (bounding)
        {
            case IEB_ERROR:
            {
                FatalErrorInFunction
                    << "value (" << lookupValue << ") greater than upper "
                    << "bound (" << maxLimit << ")\n"
                    << exit(FatalError);
                break;
            }
            case IEB_WARN:
            {
                WarningInFunction
                    << "value (" << lookupValue << ") greater than upper "
                    << "bound (" << maxLimit << ")\n"
                    << "    Continuing with the last entry" << endl;

                // Behaviour as per 'CLAMP'
                return list.last().second();
            }
            case IEB_CLAMP:
            {
                return list.last().second();
            }
            case IEB_REPEAT:
            {
                // Adjust lookupValue <= maxLimit
                const scalar span = maxLimit - minLimit;
                lookupValue = fmod(lookupValue - minLimit, span) + minLimit;
                break;
            }
            case IEB_EXTRAPOLATE:
            {
                // Linear extrapolation using last two points
                const scalar xNm1 = list[n-2].first();
                const scalar xN   = list[n-1].first();
                const Type&  yNm1 = list[n-2].second();
                const Type&  yN   = list[n-1].second();

                const scalar dx  = xN - xNm1;
                const scalar dxi = 1.0/dx;

                return
                (
                    yN
                  + (yN - yNm1)*(lookupValue - xN)*dxi
                );
            }
        }
    }

    label lo = 0;
    label hi = 0;

    // Look for the correct range
    for (label i = 0; i < n; ++i)
    {
        if (lookupValue >= list[i].first())
        {
            lo = hi = i;
        }
        else
        {
            hi = i;
            break;
        }
    }

    if (lo == hi)
    {
        return list[hi].second();
    }
    else if (hi == 0)
    {
        // This treatment should only occur under these conditions:
        //  -> the 'REPEAT' treatment
        //  -> (0 <= value <= minLimit)
        //  -> minLimit > 0
        // Use the value at maxLimit as the value for value=0
        lo = n - 1;

        return
        (
            list[lo].second()
          + (list[hi].second() - list[lo].second())
          * (lookupValue / minLimit)
        );
    }

    // Normal interpolation
    return
    (
        list[lo].second()
      + (list[hi].second() - list[lo].second())
      * (lookupValue - list[lo].first())
      / (list[hi].first() - list[lo].first())
    );
}


template<class Type>
Type interExtrapolationTable<Type>::interpolateValue
(
    scalar lookupValue
) const
{
    return interpolateValue(*this, lookupValue, bounding_);
}


template<class Type>
tmp<Field<Type>>
interExtrapolationTable<Type>::interpolateValues
(
    const UList<scalar>& vals
) const
{
    auto tfld = tmp<Field<Type>>::New(vals.size());
    auto& fld = tfld.ref();

    forAll(fld, i)
    {
        fld[i] = interpolateValue(vals[i]);
    }

    return tfld;
}


// * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * * * * //

template<class Type>
void interExtrapolationTable<Type>::operator=
(
    const interExtrapolationTable<Type>& rhs
)
{
    if (this == &rhs)
    {
        return;
    }

    static_cast<List<value_type>&>(*this) = rhs;
    bounding_ = rhs.bounding_;
    fileName_ = rhs.fileName_;
    reader_.reset(rhs.reader_.clone());
}


template<class Type>
const Tuple2<scalar, Type>&
interExtrapolationTable<Type>::operator[](label idx) const
{
    const List<value_type>& list = *this;
    const label n = list.size();

    if (n <= 1)
    {
        idx = 0;

        #ifdef FULLDEBUG
        if (!n)
        {
            FatalErrorInFunction
                << "Cannot interpolate from zero-sized table" << nl
                << exit(FatalError);
            break;
        }
        #endif
    }
    else if (idx < 0)
    {
        switch (bounding_)
        {
            case IEB_ERROR:
            {
                FatalErrorInFunction
                    << "index (" << idx << ") underflow" << nl
                    << exit(FatalError);
                break;
            }
            case IEB_WARN:
            {
                WarningInFunction
                    << "index (" << idx << ") underflow" << nl
                    << "    Continuing with the first entry" << nl;

                // Behaviour as per 'CLAMP'
                idx = 0;
                break;
            }
            case IEB_CLAMP:
            {
                idx = 0;
                break;
            }
            case IEB_REPEAT:
            {
                while (idx < 0)
                {
                    idx += n;
                }
                break;
            }
            case IEB_EXTRAPOLATE:
            {
                // EXTRAPOLATE does not make sense for indices;
                // treat like CLAMP to first entry
                idx = 0;
                break;
            }
        }
    }
    else if (idx >= n)
    {
        switch (bounding_)
        {
            case IEB_ERROR:
            {
                FatalErrorInFunction
                    << "index (" << idx << ") overflow" << nl
                    << exit(FatalError);
                break;
            }
            case IEB_WARN:
            {
                WarningInFunction
                    << "index (" << idx << ") overflow" << nl
                    << "    Continuing with the last entry" << nl;

                // Behaviour as per 'CLAMP'
                idx = n - 1;
                break;
            }
            case IEB_CLAMP:
            {
                idx = n - 1;
                break;
            }
            case IEB_REPEAT:
            {
                while (idx >= n)
                {
                    idx -= n;
                }
                break;
            }
            case IEB_EXTRAPOLATE:
            {
                // EXTRAPOLATE does not make sense for indices;
                // treat like CLAMP to last entry
                idx = n - 1;
                break;
            }
        }
    }

    return list[idx];
}


template<class Type>
Type interExtrapolationTable<Type>::operator()(scalar lookupValue) const
{
    return interpolateValue(*this, lookupValue, bounding_);
}


// -------------------------------------------------------------------------
// Explicit template instantiation
// -------------------------------------------------------------------------

#include "scalar.H"
#include "vector.H"

template class interExtrapolationTable<scalar>;
template class interExtrapolationTable<vector>;


} // End namespace Foam

// ************************************************************************* //
