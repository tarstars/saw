/// @file errors.h
///
/// @author Dmitry Azhichakov <dmitry@dsa.pp.ru>

#ifndef __ERRORS_H__INCLUDED__
#define __ERRORS_H__INCLUDED__

#include <exception>


namespace farn {

    class Error: public std::exception {
    public:
        Error(const char* what_msg) :
            what_(what_msg)
        {}

        virtual const char* what() const throw() {
            return what_;
        }
    private:
        const char* what_;
    };

}

#endif // __ERRORS_H__INCLUDED__
