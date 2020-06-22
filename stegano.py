import numpy as np
from PIL import Image
import argparse


def read_image_greyscale(image_file):
    """
    Translates greyscale png image to np.array vector of bits

    :param image_file: image (*png)
    :return: tuple of shape, vector of bits describing the image
    """
    data = np.array(Image.open(image_file).convert('L'))
    shape = data.shape
    data = data.flatten()
    return shape, np.unpackbits(data).astype('bool_')


def write_image_greyscale(x, y, bit_vector, output_image_file):
    """
    Translates vector of bits to greyscale png image sized (x,y)

    :param x: np.shape[0]
    :param y: np.shape[1]
    :param bit_vector: np.array dtype='bool_' vector
    :param output_image_file: location for created image (*png)
    """
    data = np.packbits(bit_vector)
    data = np.reshape(data, (x, y))
    Image.fromarray(data).save(output_image_file)


def read_DNA(text_file):
    """
    Translates fasta file to vector of bits with coding:
    A=00
    C=01
    G=10
    T=11

    :param text_file: fasta file (*txt)
    :return: vector of bits

    """
    with open(text_file, 'r') as f:
        data = ''.join(f.readlines()[1:]).replace('\n', '')  # getting DNA code from text_file (skipping the first line)
    vec = np.zeros((2*len(data)), dtype='bool_')  # every letter in encoding has 2 bits
    for i in range(len(data)):
        if data[i] == 'C':
            vec[2*i+1] = 1
        elif data[i] == 'G':
            vec[2*i] = 1
        elif data[i] == 'T':
            vec[2*i] = 1
            vec[2*i+1] = 1
    return vec


def write_DNA(bit_vector, output_text_file):
    """
    Translates np.array vector of bits to txt file of DNA with coding:
    A=00
    C=01
    G=10
    T=11

    :param bit_vector: np.array dtype = 'bool_' vector
    :param output_text_file: location for created text file (*txt)
    """
    with open(output_text_file, 'w+') as f:
        for i in range(int(len(bit_vector)/2)):  # 2 bits = 1 character
            if not bit_vector[2*i] and not bit_vector[2*i+1]:
                f.write('A')
            elif not bit_vector[2*i] and bit_vector[2*i+1]:
                f.write('C')
            elif bit_vector[2*i] and not bit_vector[2*i+1]:
                f.write('G')
            elif bit_vector[2*i] and bit_vector[2*i+1]:
                f.write('T')


def encode_DNA(bit_vector, image_file, output_file):
    """
    Encode vector of bits representing DNA sequence into last bits of bytes in an array representation
    (dtype='uint8') of png image. First coded bit describe whether encoded vector represents DNA sequence (0) or png
    image (1). Next 32 bits (in DNA encoding) are bit representation of uint32 length of DNA code.

    :param bit_vector: np.array dtype = 'bool_' vector
    :param image_file: image in witch one want to encode the DNA (*png)
    :param output_file: location for created image (*png)
    :return:
    """
    im = np.array(Image.open(image_file))
    shape = im.shape
    if shape[0]*shape[1] < len(bit_vector)+33:  # checking if image is big enough for encoding
        raise ValueError('image file is too small to encode the message')
    header = np.zeros((33), dtype = 'bool_')
    header[1:] = np.array([int(x) for x in '{0:032b}'.format(int(len(bit_vector)/2))])  # 2 bits = 1 character
    im = im.reshape((shape[0]*shape[1]*3))  # flattening the image
    for i in range(len(header)+ len(bit_vector)):
        if i < 33:  # encoding the header (type and length)
            if im[i]%2 != header[i]:
                if im[i] == 0:
                    im[i] += 1
                else:
                    im[i] -= 1
        else:  # encoding the DNA sequence
            if im[i]%2 != bit_vector[i-33]:
                if im[i] == 0:
                    im[i] += 1
                else:
                    im[i] -= 1
    Image.fromarray(im.reshape(shape)).save(output_file)


def encode_image(bit_vector, x, y, image_file, output_file):
    """
    Encode vector of bits representing a greyscale image sized (x,y) into last bits of bytes in an array representation
    (dtype='uint8') of png image. First coded bit describe whether encoded vector represents DNA sequence (0) or png
    image (1). Next 32 bits (in png image encoding) are uint16 bit representations of respectively x and y (image size).
    :param bit_vector: np.array dtype = 'bool_' vector
    :param x: np.shape[0]
    :param y: np.shape[1]
    :param image_file: image in witch one want to encode the DNA (*png)
    :param output_file: location for created image (*png)
    :return:
    """
    im = np.array(Image.open(image_file))
    shape = im.shape
    if shape[0]*shape[1] < len(bit_vector)+33:  # checking if image is big enough for encoding
        raise ValueError('image file is too small to encode the message')
    header = np.ones((33), dtype = 'bool_')
    header[1:17] = np.array([int(x) for x in '{0:016b}'.format(x)])
    header[17:] = np.array([int(x) for x in '{0:016b}'.format(y)])
    im = im.reshape((shape[0]*shape[1]*3))
    for i in range(len(header)+ len(bit_vector)):
        if i < 33:  # encoding the header
            if im[i]%2 != header[i]:
                if im[i] == 0:
                    im[i] += 1
                else:
                    im[i] -= 1
        else:    # encoding bit representation of greyscale image
            if im[i]%2 != bit_vector[i-33]:
                if im[i] == 0:
                    im[i] += 1
                else:
                    im[i] -= 1
    Image.fromarray(im.reshape(shape)).save(output_file)


def get_header(flat_image_array):
    """
    Gets first 33 encoded bits (header)
    :param flat_image_array: np.array vector
    :return: np.array vector header
    """
    header = np.ones((33), dtype='bool_')
    for i in range(33):
        if flat_image_array[i]%2 == 0:
            header[i] = 0
    return header


def get_bit_vector(flat_image_array, length):
    """
    Gest encoded vector of bits
    :param flat_image_array: np.array vector
    :param length: length of mined vector of bits
    """
    bit_vector = np.ones((length), dtype='bool_')
    for i in range(length):
        if flat_image_array[i+33]%2 == 0:  # skipping the header
            bit_vector[i] = 0
    return bit_vector


def decode(image_file, output_file):
    """
    Decodes hidden image file or DNA sequence and saves it to the output_file location (image -> (*png), DNA -> (*txt))
    :param image_file: location of encoded image (*png)
    :param output_file: location for decoded file (*png) or (*txt)
    """
    im = np.array(Image.open(image_file)).flatten()
    header = get_header(im)
    if header[0] == 0:  # encoded file is a DNA sequence
        b32 = np.flip(header[1:])
        length = 0
        for i in range(len(b32)):  # encoding length of the DNA sequence
            length += 2**i*b32[i]
        write_DNA(get_bit_vector(im,length*2), output_file)  # encoding and saving the sequence (2 bits = 1 character)
    elif header[0] == 1:  # encoded file is a greyscale image
        b16x = np.flip(header[1:17])
        b16y = np.flip(header[17:])
        x, y = 0, 0
        for i in range(len(b16x)):  # encoding the size of the image
            x += 2**i*b16x[i]
            y += 2**i*b16y[i]
        write_image_greyscale(x, y, get_bit_vector(im,x*y*8), output_file)  # encoding and saving the greyscale image (8 bits = 1 byte)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("image_in", type=str,
                        help="base image for encoding/decoding")
    parser.add_argument("file_out", type=str,
                        help="location for saving encoded/decoded file")
    parser.add_argument("--ed", "--encodeDNA", type=str, default="",
                        help="encodes provided DNA sequence (*txt) into an image")
    parser.add_argument("--ei", "--encodeimage", type=str, default="",
                        help="encodes provided greyscale image (*png) into an image")
    parser.add_argument("--d", "--decode", action="store_true",
                        help="decodes provided image_in and saves decoded file in file_out")
    args = parser.parse_args()
    if args.ed != "":
        encode_DNA(read_DNA(args.ed), args.image_in, args.file_out)
    elif args.ei != "":
        res = read_image_greyscale(args.ei)
        encode_image(res[1], res[0][0], res[0][1], args.image_in, args.file_out)
    elif args.d:
        decode(args.image_in, args.file_out)


main()