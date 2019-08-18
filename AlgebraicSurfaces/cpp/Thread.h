//
//  Thread.h
//  HillClimbImage
//
//  Created by asd on 29/04/2019.
//  Copyright © 2019 voicesync. All rights reserved.
//

#ifndef Thread_h
#define Thread_h

#include <thread>

using std::thread;

class Thread {
public:
    Thread(int size) : nth(thread::hardware_concurrency()), segSz(size/nth),
     threads(new thread[nth]), size(size) {}
    
    ~Thread() {
        delete[]threads;
    }
    static int getnthreads() {return thread::hardware_concurrency(); }
    
    int from(int t) { return t*segSz; }
    int to(int t) { return ((t==nth-1) ? size : (t+1)*segSz); }
    
    void run(std::function<void(int, int, int)> const& lambda) {
        for (int t=0; t<nth; t++) {
            threads[t]=thread([this, lambda, t](){
                lambda(t, from(t), to(t));
            });
        }
        for (int t=0; t<nth; t++) threads[t].join();
    }
    void run(std::function<void(int)> const& lambda) {
        for (int t=0; t<nth; t++) {
            threads[t]=thread([this, lambda, t](){
                for (int i=from(t); i<to(t); i++)
                    lambda(i);
            });
        }
        for (int t=0; t<nth; t++) threads[t].join();
    }
    void run(std::function<void(void)> const& lambda) {
        for (int t=0; t<nth; t++) {
            threads[t]=thread([this, lambda, t](){
                for (int i=from(t); i<to(t); i++)
                    lambda();
            });
        }
        for (int t=0; t<nth; t++) threads[t].join();
    }
    int nth, segSz;
    thread *threads;
    int size;
    
};

#endif /* Thread_h */
