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

\*---------------------------------------------------------------------------*/

#include "mappedLeafTempFvPatchScalarField.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "indexedOctree.H"
#include "treeDataPoint.H"
//v8: #include "mappedFvPatchBaseBase.H"
#include "mappedInternalPatchBase.H" //v12: mappedInternal uses separate class hierarchy

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::mappedLeafTempFvPatchScalarField::
mappedLeafTempFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(p, iF),
    fieldName_(iF.name())
    //v8: mapperPtr_(nullptr)
{}


Foam::mappedLeafTempFvPatchScalarField::
mappedLeafTempFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchScalarField(p, iF, dict),
    fieldName_(dict.lookupOrDefault<word>("field", iF.name()))
{
    if (!isA<mappedInternalPatchBase>(p.patch()))
    {
        FatalIOErrorInFunction(dict)
            << "Field " << iF.name() << " on patch " << p.name()
            << " is not of mappedInternal type"
            << exit(FatalIOError);
    }
}


Foam::mappedLeafTempFvPatchScalarField::
mappedLeafTempFvPatchScalarField
(
    const mappedLeafTempFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchScalarField(ptf, p, iF, mapper),
    fieldName_(ptf.fieldName_)
{}


Foam::mappedLeafTempFvPatchScalarField::
mappedLeafTempFvPatchScalarField
(
    const mappedLeafTempFvPatchScalarField& ptf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    fixedValueFvPatchScalarField(ptf, iF),
    fieldName_(ptf.fieldName_)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

const Foam::mappedInternalPatchBase&
Foam::mappedLeafTempFvPatchScalarField::mapper() const
{
    return refCast<const mappedInternalPatchBase>(this->patch().patch());
}


void Foam::mappedLeafTempFvPatchScalarField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    const int oldTag = UPstream::msgType();
    UPstream::msgType() = oldTag + 1;

    const mappedInternalPatchBase& mipb = this->mapper();
    const fvMesh& airMesh = refCast<const fvMesh>(mipb.nbrMesh());
    const volScalarField& Tl = airMesh.lookupObject<volScalarField>("Tl");
    const volScalarField& LAD = airMesh.lookupObject<volScalarField>("LAD");

    scalarField& Tp = *this;

    List<List<point>> vegCellCentres(Pstream::nProcs());
    List<List<scalar>> vegCellValues(Pstream::nProcs());
    const scalarField& TlI = Tl.primitiveField();
    const scalarField& LADI = LAD.primitiveField();
    const vectorField& CI = airMesh.cellCentres();
    forAll(TlI, cellI)
    {
        if (LADI[cellI] > 0 && TlI[cellI] > 0)
        {
            vegCellCentres[Pstream::myProcNo()].append(CI[cellI]);
            vegCellValues[Pstream::myProcNo()].append(TlI[cellI]);
        }
    }
    Pstream::gatherList(vegCellCentres);
    Pstream::scatterList(vegCellCentres);

    Pstream::gatherList(vegCellValues);
    Pstream::scatterList(vegCellValues);

    //flatten canopy cells of all processors
    DynamicList<point> allCentres;
    DynamicList<scalar> allValues;
    forAll(vegCellCentres, proci)
    {
        allCentres.append(vegCellCentres[proci]);
        allValues.append(vegCellValues[proci]);
    }

    //nearest canopy cell via octree (keep the current value when no
    //LAD>0/Tl>0 cells exist globally yet, e.g. before the first vegetation solve)
    if (allCentres.size())
    {
        pointField centres;
        centres.transfer(allCentres);

        //fixed padding so a single canopy cell does not give a zero-size box
        treeBoundBox bb(centres);
        bb.min() -= vector::uniform(1e-3);
        bb.max() += vector::uniform(1e-3);

        const indexedOctree<treeDataPoint> tree
        (
            treeDataPoint(centres),
            bb,
            8,      // maxLevel
            10,     // leafsize
            3.0     // duplicity
        );

        const vectorField& Cf = this->patch().Cf();
        forAll(Tp, i)
        {
            const pointIndexHit hit = tree.findNearest(Cf[i], sqr(great));
            if (hit.hit())
            {
                Tp[i] = allValues[hit.index()];
            }
        }
    }

    UPstream::msgType() = oldTag;

    fixedValueFvPatchScalarField::updateCoeffs();
}


void Foam::mappedLeafTempFvPatchScalarField::write(Ostream& os) const
{
    fvPatchScalarField::write(os);
    writeEntry(os, "field", fieldName_);
    writeEntry(os, "value", *this);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
    makePatchTypeField
    (
        fvPatchScalarField,
        mappedLeafTempFvPatchScalarField
    );
}

// ************************************************************************* //
