#ifndef ESCAPE_ERRORCODE_H_
#define ESCAPE_ERRORCODE_H_

namespace Escape
{
  enum ErrorCode : int
  {
    ecNone
    , ecInvalidInput
    , ecSystemError
    , ecUnsupportedFormat
    , ecIOError
  };
}
#endif
