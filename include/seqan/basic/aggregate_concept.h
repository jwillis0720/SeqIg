// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_AGGREGATE_CONCEPT_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_AGGREGATE_CONCEPT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @concept AggregateConcept
 * 
 * @brief Aggregate types contain a fixed number of fixed-size values (pairs, triples, tuples).
 * 
 * @section Remarks
 * 
 * Stream output operators are not shown in the function list below, but required.
 * 
 * Comparison operators are not shown in the function list below, but required.
 */

/*!
 * @fn AggregateConcept#operator<<
 * @brief Stream output operator.
 *
 * @signature TStream AggregateConcept::operator<<(stream, aggregate);
 *
 * @param stream    The <tt>std::ostream</tt> to write to.
 * @param aggregate The aggregate type to write to the stream.
 *
 * @return TStream Reference to <tt>stream</tt> after writing <tt>aggregate</tt> to it.
 */

/**
.Concept.AggregateConcept
..summary:Aggregate types contain a fixed number of fixed-size values.
..remarks:Stream output operators are not shown in the function list below, but required.
..remarks:Comparison operators are not shown in the function list below, but required.

.Function.clear.concept:Concept.AggregateConcept
.Function.value.concept:Concept.AggregateConcept
.Function.assignValue.concept:Concept.AggregateConcept

.Metafunction.LENGTH.concept:Concept.AggregateConcept
.Metafunction.Value.concept:Concept.AggregateConcept
 */

/*!
 * @defgroup AggregateTags Aggregate Tags
 * @brief Tags to use in aggregate (e.g. Pair, Triple, and Tuple) types.
 */

/*!
 * @tag AggregateTags#Pack
 * @headerfile <seqan/basic.h>
 * @brief Tag to mark a packed specialization that disables address alignment for members.
 *
 * @signature typedef Tag<Pack_> Pack;
 */

/**
.Tag.Pack
..cat:Aggregates
..summary:Tag to mark a packed specialization that disables address alignment for members.
..signature:Pack
..include:seqan/basic.h
 */

struct Pack_;
typedef Tag<Pack_> Pack;

// TODO(holtgrew): We need @tparam for tag in the Dox system.

/*!
 * @tag AggregateTags#BitPacked
 * @headerfile <seqan/basic.h>
 * @brief Tag to mark a bit-packed specialization that avoids to waste bits.
 *
 * @signature template <[unsinged BITSIZE1[, unsigned BITSIZE2]]>
 *            struct BitPacked;
 *
 * BITSIZE1 The number of bits for the first entry.
 *
 * BITSIZE2 The number of bits for the second entry.
 */

/**
.Tag.BitPacked
..cat:Aggregates
..summary:Tag to mark a bit-packed specialization that avoids to waste bits.
..signature:BitPacked<BITSIZE1, BITSIZE2>
..param.BITSIZE1:Number of bits used for first element.
...type:nolink:$unsigned$
..param.BITSIZE2:Number of bits used for second element.
...type:nolink:$unsigned$
..include:seqan/basic.h
 */

template <unsigned BITSIZE1 = 16, unsigned BITSIZE2 = 16>
struct BitPacked;

// ============================================================================
// Metafunctions
// ============================================================================

/*!
 * @mfn MakePacked
 * @headerfile <seqan/basic.h>
 * @brief Return the corresponding packed type for a type.
 *
 * @signature MakePacked<TAggregate>::Type;
 *
 * @tparam TAggregate The aggregate type to transform.
 *
 * @return Type The resulting packed type.
 */

/**
.Metafunction.MakePacked
..cat:Aggregates
..summary:Return the corresponding packed type of a type.
..signature:MakePacked<TAggregate>
..param.TAggregate:An aggregate type.
..returns:The corresponding packed aggregate.
..include:seqan/basic.h
 */

template <typename T>
struct MakePacked
{
    typedef T Type;
};

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_AGGREGATE_CONCEPT_H_
