/*
 * Replicated.cpp
 *
 */

#ifndef PROTOCOLS_MULTIREPLICATED_HPP_
#define PROTOCOLS_MULTIREPLICATED_HPP_

#include "MultiReplicated.h"
#include "Processor/Processor.h"
#include "Tools/benchmarking.h"

#include "SemiShare.h"
#include "SemiMC.h"
#include "ReplicatedInput.h"
#include "Rep3Share2k.h"

#include "SemiMC.hpp"
#include "Math/Z2k.hpp"


template<class T>
MultiReplicated<T>::MultiReplicated(Player& P) : MultiReplicatedBase(P)
{
    assert(T::length == 2);
}

template<class T>
MultiReplicated<T>::MultiReplicated(const MultiReplicatedBase& other) :
        MultiReplicatedBase(other)
{
}

inline MultiReplicatedBase::MultiReplicatedBase(Player& P) : P(P)
{
    assert(P.num_players() == 3);
    if (not P.is_encrypted())
        insecure("unencrypted communication");
    if (P.my_num() == 0){
        for (int i = 0; i < 2; ++i) {
            shared_prngs[i].ReSeed();
            octetStream os;
            os.append(shared_prngs[i].get_seed(), SEED_SIZE);
            P.send_relative(1 - 2*i, os);
        }
    } else {
        octetStream os;
        P.receive_player(0, os, true);
        shared_prngs[0].SetSeed(os.get_data());
    }
}

inline MultiReplicatedBase::MultiReplicatedBase(Player& P, array<PRNG, 2>& prngs) :
        P(P)
{
    for (int i = 0; i < 2; i++)
        shared_prngs[i].SetSeed(prngs[i]);
}

inline MultiReplicatedBase MultiReplicatedBase::branch()
{
    return {P, shared_prngs};
}

template<class T>
void MultiReplicated<T>::init_mul(SubProcessor<T>* proc)
{
    (void) proc;
    init_mul();
}

template<class T>
void MultiReplicated<T>::init_mul(Preprocessing<T>& prep, typename T::MAC_Check& MC)
{
    (void) prep, (void) MC;
    init_mul();
}

template<class T>
void MultiReplicated<T>::init_mul()
{
    if(P.my_num() == 0)
        os.resize(1);
    else if (P.my_num() == 1)
        os.resize(2);
    else
        os.resize(3);
    for (auto& o : os)
        o.reset_write_head();
    add_shares.clear();
}


template<class T>
inline typename T::clear MultiReplicated<T>::prepare_mul(const vector<T>& operands, int n)
{

    T x = operands[0];
    T y = operands[1];

    typename T::value_type comm_elem;
    if (P.my_num() == 0){
        comm_elem = x.local_mul_0(y);
        randomize_comm_elem_0(comm_elem, n);
    } else if (P.my_num() == 1){
        comm_elem = x.local_mul_1(y);
        randomize_comm_elem_1(comm_elem, n);
    } else {
        comm_elem = x.local_mul_2(y);
        randomize_comm_elem_2(comm_elem, n);
    }

    return comm_elem;
}

template<class T>
inline void MultiReplicated<T>::randomize_comm_elem_0(const typename T::clear& comm_elem,
                                           int n)
{
    auto elem = comm_elem;
    typename T::clear rnd[3];
    for (int i = 0; i < 2; i++){
        rnd[i].randomize(shared_prngs[0], n);
        next_relative_rnds.push_back(rnd[i]);
    }
    rnd[2].randomize(shared_prngs[1], n);
    prev_relative_rnds.push_back(rnd[2]);
    add_shares.push_back(elem);
//    cout << "party 1 saved " << elem << endl;
    elem += rnd[1];
//    cout << "party 1 sent " << elem << endl;
    elem.pack(os[0], n);
}

template<class T>
inline void MultiReplicated<T>::randomize_comm_elem_1(const typename T::clear& comm_elem,
                                                      int n)
{
    auto elem = comm_elem;
    typename T::clear rnd[2];
    for (int i = 0; i < 2; i++){
        rnd[i].randomize(shared_prngs[0], n);
        prev_relative_rnds.push_back(rnd[i]);
    }
    add_shares.push_back(elem);
    elem += rnd[0];
//    cout << "party 2 sent " << elem << endl;
    elem.pack(os[0], n);
}

template<class T>
inline void MultiReplicated<T>::randomize_comm_elem_2(const typename T::clear& comm_elem,
                                                      int n)
{
    auto elem = comm_elem;
    typename T::clear rnd;
    rnd.randomize(shared_prngs[0], n);
    next_relative_rnds.push_back(rnd);
    add_shares.push_back(elem);
    elem += rnd;
//    cout << "party 3 sent " << elem << endl;
    elem.pack(os[0], n);
}

template<class T>
void MultiReplicated<T>::exchange()
{
    if (P.my_num() == 0)
        P.send_relative(-1, os[0]);
    if (P.my_num() == 1){
        P.send_relative(1, os[0]);
        P.receive_relative(1, os[1]);
    }
    if (P.my_num() == 2){
        P.send_relative(-1, os[0]);
        P.receive_relative(1, os[1]);
        P.receive_relative(-1, os[2]);
    }
}

template<class T>
void MultiReplicated<T>::start_exchange()
{
    P.send_relative(1, os[0]);
}

template<class T>
void MultiReplicated<T>::stop_exchange()
{
    P.receive_relative(-1, os[1]);
}

template<class T>
inline T MultiReplicated<T>::finalize_mul(int n)
{
    T result;
    if (P.my_num() == 0){
        vector<typename T::clear> triple_rnd = {prev_relative_rnds.next(), next_relative_rnds.next(),
                                                next_relative_rnds.next()};
//        cout << "alpha " << triple_rnd[0] << " beta " << triple_rnd[1] << " gamma " << triple_rnd[2] << endl;
        typename T::clear tmp = add_shares.next();
//        cout << "party 1 should have saved " << tmp << endl;
        result[0] = tmp + triple_rnd[0] + triple_rnd[2];
        result[1] = triple_rnd[1] + triple_rnd[2];
//        cout << "party 1 shares " << result[0] << " " << result[1] << endl;
    }
    else if (P.my_num() == 1){
        typename T::clear received;
        vector<typename T::clear> pair_rnd = {prev_relative_rnds.next(), prev_relative_rnds.next()};
        received.unpack(os[1], n);
//        cout << " beta " << pair_rnd[0] << " gamma " << pair_rnd[1] << endl;
//        cout << "party 2 received " << received << endl;
        result[0] = received + add_shares.next() + pair_rnd[1];
        result[1] = pair_rnd[0] + pair_rnd[1];
//        cout << "party 2 shares " << result[0] << " " << result[1] << endl;
    }
    else if (P.my_num() == 2){
        typename T::clear received_prev, received_next;
        typename T::clear rnd = next_relative_rnds.next();
        received_next.unpack(os[1], n);
        received_prev.unpack(os[2], n);
//        cout << " alpha " << rnd << endl;
//        cout << "party 3 received from 1 " << received_next << " and from 2 " << received_prev << endl;
        result[0] = add_shares.next() + received_next + received_prev;
        result[1] = received_next + rnd;
//        cout << "party 3 shares " << result[0] << " " << result[1] << endl;
    }

    return result;
}

template<class T>
inline void MultiReplicated<T>::init_dotprod(SubProcessor<T>* proc)
{
    init_mul(proc);
    dotprod_share.assign_zero();
}

template<class T>
inline void MultiReplicated<T>::prepare_dotprod(const T& x, const T& y)
{
    dotprod_share += x.local_mul(y);
}

template<class T>
inline void MultiReplicated<T>::next_dotprod()
{
    prepare_reshare(dotprod_share);
    dotprod_share.assign_zero();
}

template<class T>
inline T MultiReplicated<T>::finalize_dotprod(int length)
{
    (void) length;
    this->counter++;
    return finalize_mul();
}

template<class T>
T MultiReplicated<T>::get_random()
{
    T res;
    for (int i = 0; i < 2; i++)
        res[i].randomize(shared_prngs[i]);
    return res;
}

template<class T>
void MultiReplicated<T>::trunc_pr(const vector<int>& regs, int size,
                             SubProcessor<T>& proc)
{
    ::trunc_pr(regs, size, proc);
}

#endif
