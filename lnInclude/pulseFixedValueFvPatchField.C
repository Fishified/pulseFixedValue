/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2019 OpenFOAM Foundation
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

#include "pulseFixedValueFvPatchField.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"

// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
Foam::scalar Foam::pulseFixedValueFvPatchField<Type>::pulseFraction() const
{
    const scalar t = this->db().time().timeOutputValue();
    
    // t = 42; period = 5; duration = 0.5
    // t/period = 42/5 = 8.4 -> 0.4 -> < 0.5
    
    scalar cycleFraction = fmod(t/period_, 1.0);
    
    if (cycleFraction > duration_)
    {
        return 0.0;
    }
    else
    {
        return 1.0;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::pulseFixedValueFvPatchField<Type>::
pulseFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(p, iF),
    period_(0.0),
    duration_(0.0),
    baseValue_(Zero),
    pulseValue_(p.size(), Zero)
{
}


template<class Type>
Foam::pulseFixedValueFvPatchField<Type>::
pulseFixedValueFvPatchField
(
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<Type>(p, iF),
    period_(readScalar(dict.lookup("period"))),
    duration_(readScalar(dict.lookup("duration"))),
    baseValue_(dict.lookupOrDefault<Type>("baseValue", Zero)),
    pulseValue_("pulseValue", dict, p.size())
{
    if (duration_ < 0 || duration_ > 1)
    {
        FatalErrorInFunction
            << "duration = " << duration_ << " is invalid (0 < duration < 1)"
            << exit(FatalError);
    }

    fixedValueFvPatchField<Type>::evaluate();

    /*
    // Initialise with the value entry if evaluation is not possible
    fvPatchField<Type>::operator=
    (
        Field<Type>("value", dict, p.size())
    );
    */
}


template<class Type>
Foam::pulseFixedValueFvPatchField<Type>::
pulseFixedValueFvPatchField
(
    const pulseFixedValueFvPatchField<Type>& ptf,
    const fvPatch& p,
    const DimensionedField<Type, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<Type>(ptf, p, iF, mapper),
    period_(ptf.period_),
    duration_(ptf.duration_),
    baseValue_(ptf.baseValue_),
    pulseValue_(ptf.pulseValue_, mapper)
{}


template<class Type>
Foam::pulseFixedValueFvPatchField<Type>::
pulseFixedValueFvPatchField
(
    const pulseFixedValueFvPatchField<Type>& ptf
)
:
    fixedValueFvPatchField<Type>(ptf),
    period_(ptf.period_),
    duration_(ptf.duration_),
    baseValue_(ptf.baseValue_),
    pulseValue_(ptf.pulseValue_)
{}


template<class Type>
Foam::pulseFixedValueFvPatchField<Type>::
pulseFixedValueFvPatchField
(
    const pulseFixedValueFvPatchField<Type>& ptf,
    const DimensionedField<Type, volMesh>& iF
)
:
    fixedValueFvPatchField<Type>(ptf, iF),
    period_(ptf.period_),
    duration_(ptf.duration_),
    baseValue_(ptf.baseValue_),
    pulseValue_(ptf.pulseValue_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::pulseFixedValueFvPatchField<Type>::autoMap
(
    const fvPatchFieldMapper& m
)
{
    fixedValueFvPatchField<Type>::autoMap(m);
    pulseValue_.autoMap(m);
}


template<class Type>
void Foam::pulseFixedValueFvPatchField<Type>::rmap
(
    const fvPatchField<Type>& ptf,
    const labelList& addr
)
{
    fixedValueFvPatchField<Type>::rmap(ptf, addr);

    const pulseFixedValueFvPatchField<Type>& tiptf =
        refCast<const pulseFixedValueFvPatchField<Type>>(ptf);

    pulseValue_.rmap(tiptf.pulseValue_, addr);
}


template<class Type>
void Foam::pulseFixedValueFvPatchField<Type>::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    fixedValueFvPatchField<Type>::operator==
    (
        baseValue_*(1.0 - pulseFraction())
      + pulseValue_*pulseFraction()
    );

    fixedValueFvPatchField<Type>::updateCoeffs();
}


template<class Type>
void Foam::pulseFixedValueFvPatchField<Type>::write
(
    Ostream& os
) const
{
    fvPatchField<Type>::write(os);
    os.writeKeyword("period") << period_ << token::END_STATEMENT << nl;
    os.writeKeyword("duration") << duration_ << token::END_STATEMENT << nl;
    os.writeKeyword("baseValue") << baseValue_ << token::END_STATEMENT << nl;
    pulseValue_.writeEntry("pulseValue", os);
    this->writeEntry("value", os);
}



// ************************************************************************* //
