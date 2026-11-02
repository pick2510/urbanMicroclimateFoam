/*---------------------------------------------------------------------------*\
Global Function
    Foam::nearestPatchField

Description
    Values of nbrFld (defined at face centres nbrCf, any region) at the
    nearest nbrCf of each point of Cf, searched over all processors.

    Replaces ad-hoc mappedPatchBase mappers: the OpenFOAM tree search aborts
    on a 1-face processor patch in older v12 builds, and patchToPatch
    "nearest" leaves faces without a nearby neighbour uncoupled (NaN).
\*---------------------------------------------------------------------------*/

#ifndef nearestPatchField_H
#define nearestPatchField_H

#include "pointField.H"
#include "scalarField.H"
#include "DynamicList.H"
#include "Pstream.H"
#include "indexedOctree.H"
#include "treeDataPoint.H"

namespace Foam
{

inline tmp<scalarField> nearestPatchField
(
    const pointField& Cf,
    const pointField& nbrCf,
    const scalarField& nbrFld
)
{
    List<List<point>> procCf(Pstream::nProcs());
    List<List<scalar>> procFld(Pstream::nProcs());
    procCf[Pstream::myProcNo()] = nbrCf;
    procFld[Pstream::myProcNo()] = nbrFld;

    Pstream::gatherList(procCf);
    Pstream::scatterList(procCf);
    Pstream::gatherList(procFld);
    Pstream::scatterList(procFld);

    DynamicList<point> allCf;
    DynamicList<scalar> allFld;
    forAll(procCf, proci)
    {
        allCf.append(procCf[proci]);
        allFld.append(procFld[proci]);
    }

    tmp<scalarField> tfld(new scalarField(Cf.size(), 0.0));
    scalarField& fld = tfld.ref();

    if (allCf.size())
    {
        pointField centres;
        centres.transfer(allCf);

        //fixed padding so a single face does not give a zero-size box
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

        forAll(fld, i)
        {
            const pointIndexHit hit = tree.findNearest(Cf[i], sqr(great));
            if (hit.hit())
            {
                fld[i] = allFld[hit.index()];
            }
        }
    }

    return tfld;
}

} // End namespace Foam

#endif

// ************************************************************************* //
