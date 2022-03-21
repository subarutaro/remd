#ifndef MACRO_H
#define MACRO_H

#ifdef __cpp_if_constexpr
#define IF_CONSTEXPR if constexpr
#else
#define IF_CONSTEXPR if
#endif

#endif
