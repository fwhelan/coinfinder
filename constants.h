//
// Created by Martin Rusilowicz on 22/08/2017.
//

#ifndef COINFINDER_ENUMS_H
#define COINFINDER_ENUMS_H

enum class ECorrection
{
        _INVALID,
        NONE,
        BONFERRONI,
        FRACTION,
};

enum class EHypothesis
{
        _INVALID,
        TWOTAILED,
        LESS,
        GREATER,
};

enum class EMethod
{
        _INVALID,
        CONNECTIVITY,
        COINCIDENCE,
};

enum class ESetMode
{
        _INVALID,
        FULL,
        INTERSECTION,
};

enum class EMaxMode
{
        _INVALID,
        ACCOMPANY,
        AVOID
};

#endif //COINFINDER_ENUMS_H
