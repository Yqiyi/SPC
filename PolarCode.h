//
// Created by Saurabh Tavildar on 5/17/16.
//

#ifndef POLARC_POLARCODE_H
#define POLARC_POLARCODE_H


#include <cstdint>
#include <vector>
#include <math.h>
#include <stack>          // std::stack

class PolarCode {


public:

    PolarCode(uint8_t num_layers, uint64_t info_length, double epsilon, uint64_t crc_size, uint64_t msg_length, std::vector<uint8_t> msg_bits) :
            _n(num_layers), _info_length(info_length), _design_epsilon(epsilon),
            _crc_size(crc_size), _msg_length(msg_length), _msg_bits(msg_bits), _llr_based_computation(true)
    {
        _block_length = (uint64_t) (1 << _n);
        _frozen_bits.resize(_block_length);
        _bit_rev_order.resize(_block_length);
        create_bit_rev_order();
        initialize_frozen_bits();
    }

    std::vector<uint8_t> encode(std::vector<uint8_t> info_bits);
	std::vector<uint8_t> encode_msg(std::vector<uint8_t> info_bits_padded);
	std::vector<uint8_t> encode_msg_ext(std::vector<uint8_t> info_bits_padded);
    std::vector<uint8_t> decode_scl_p1(std::vector<double> p1, std::vector<double> p0, uint64_t list_size);
    std::vector<uint8_t> decode_scl_llr(std::vector<double> llr, uint64_t list_size);
	//std::vector<uint8_t> spc_binary(std::vector<uint8_t> cover_bits, uint64_t cover_length, std::vector<uint8_t> msg_bits, uint64_t msg_length,  std::vector<double> prob);

    //std::vector<std::vector<double>> get_bler_quick(std::vector<double> ebno_vec, std::vector<uint8_t> list_size);

private:

    uint8_t _n;
    uint64_t _info_length;
    uint64_t _block_length;
    uint64_t _crc_size;
	uint64_t _msg_length;
	std::vector<uint8_t> _msg_bits;

    double _design_epsilon;

    std::vector<uint8_t> _frozen_bits;
    std::vector<uint64_t> _channel_order_descending;
    std::vector<std::vector<uint8_t>> _crc_matrix;
    std::vector<uint64_t> _bit_rev_order;

    void initialize_frozen_bits();
    void create_bit_rev_order();

    std::vector<uint8_t> decode_scl();
    bool _llr_based_computation;

    std::vector<std::vector<double *>> _arrayPointer_LLR;
    std::vector<double> _pathMetric_LLR;

    uint64_t _list_size;

    std::stack<uint64_t> _inactivePathIndices;
    std::vector<uint64_t > _activePath;
    std::vector<std::vector<double *>> _arrayPointer_P;
    std::vector<std::vector<uint8_t *>> _arrayPointer_C;
    std::vector<uint8_t *> _arrayPointer_Info;
    std::vector<std::vector<uint64_t>> _pathIndexToArrayIndex;
    std::vector<std::stack<uint64_t>> _inactiveArrayIndices;
    std::vector<std::vector<uint64_t>> _arrayReferenceCount;

    void initializeDataStructures();
    uint64_t assignInitialPath();
    uint64_t clonePath(uint64_t);
    void killPath(uint64_t l);

    double * getArrayPointer_P(uint64_t lambda, uint64_t  l);
    double * getArrayPointer_LLR(uint64_t lambda, uint64_t  l);
    uint8_t * getArrayPointer_C(uint64_t lambda, uint64_t  l);

    void recursivelyCalcP(uint64_t lambda, uint64_t phi);
    void recursivelyCalcLLR(uint64_t lambda, uint64_t phi);
    void recursivelyUpdateC(uint64_t lambda, uint64_t phi);

    void continuePaths_FrozenBit(uint64_t phi);
    void continuePaths_UnfrozenBit(uint64_t phi);

    uint64_t findMostProbablePath(bool check_crc);

    bool crc_check(uint8_t * info_bits_padded);

};


#endif //POLARC_POLARCODE_H
