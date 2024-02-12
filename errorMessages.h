// errorMessages.h

#ifndef ERROR_MESSAGES_H
#define ERROR_MESSAGES_H

#include "errorCode.h"
#include <string>

const std::map<ErrorCode, std::string_view> errorMessages {
        {ErrorCode::NotEnoughArguments, "You must enter the command as follows : ./render -P (or -S) <FILEPATH/file.csv>"},
        {ErrorCode::FileNotFound, "The file mentioned does not exist."},
        {ErrorCode::NotEnoughFields, "The file you mentioned does not have enough fields to plot in 3D"}
};

#endif // ERROR_MESSAGES_H