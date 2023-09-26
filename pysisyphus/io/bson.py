# [1] https://bsonspec.org/spec.html


import io
import struct


def expect(expt, handle, fmt, size):
    data = struct.unpack(fmt, handle.read(size))
    assert data == expt


def expect_byte(expt, handle):
    expect((expt,), handle, "<b", 1)


def read_int32(handle) -> int:
    """Read 4-byte integer."""
    return struct.unpack("<i", handle.read(4))[0]


def read_float64(handle) -> float:
    """Read 8-byte float/double."""
    return struct.unpack("<d", handle.read(8))[0]


def read_string(handle) -> str:
    """Read string of known length."""
    size = read_int32(handle) - 1
    string_bytes = handle.read(size)
    expect_byte(0, handle)
    return string_bytes.decode()


def read_ename(handle) -> str:
    """Read (0-terminated) string of unknown length."""
    ename_bytes = list()
    while True:
        item = handle.read(1)
        if item == b"\x00":
            break
        ename_bytes.append(item)
    return b"".join(ename_bytes).decode()


def parse_element(handle, data, size):
    kind = struct.unpack("<b", handle.read(1))[0]
    # Return early when end of object is encountered
    if kind == 0x0:
        return

    ename = read_ename(handle)

    # As of ORCA 5.0.4, only elements appearing in ORCA BSON files are implemented
    if kind == 0x1:
        value = read_float64(handle)
    elif kind == 0x2:
        value = read_string(handle)
    # Embedded document
    elif kind == 0x3:
        value = {}
        parse_document(handle, value, size)
    # Array
    elif kind == 0x4:
        value = {}
        parse_document(handle, value, size)
        # Drop integer keys and convert to list
        value = [value[str(i)] for i in range(len(value))]
    elif kind == 0x10:
        value = read_int32(handle)
    else:
        raise NotImplementedError(f"Element {kind} is not implemented!")

    data[ename] = value
    return value


def parse_document(handle, data, size=None):
    size = read_int32(handle)
    while parse_element(handle, data, size) is not None:
        pass


def load(handle) -> dict:
    """Load subset of BSON from byte-stream."""
    data = {}
    parse_document(handle, data)
    return data


def loads(bytes: bytes):
    """Load subset of BSON from bytes."""
    return load(io.BytesIO(bytes))
