// color map generator
// calcColor(u, 0, TWOPI, colourmap)

#include "color_map.h"

Color color_map(float v, float vmin, float vmax, int type) {
    float dv, vmid, r;
    Color c = {1, 1, 1};
    Color c1, c2, c3;
    float ratio;

    if (vmax < vmin) {
        dv = vmin;
        vmin = vmax;
        vmax = dv;
    }
    if (vmax - vmin < 0.000001) {
        vmin -= 1;
        vmax += 1;
    }

    if (v < vmin)     v = vmin;
    if (v > vmax)     v = vmax;
    dv = vmax - vmin;

    switch (type) {
        case 1:
            if (v < (vmin + 0.25 * dv)) {
                c = {0, 4 * (v - vmin) / dv, 1};
            }
            else if (v < (vmin + 0.5 * dv)) {
                c = {0, 1, 1 + 4 * (vmin + 0.25f * dv - v) / dv};
            }
            else if (v < (vmin + 0.75 * dv)) {
                c = { 4 * (v - vmin - 0.5f * dv) / dv, 1, 0 };
            }
            else {
                c = { 1, 1 + 4 * (vmin + 0.75f * dv - v) / dv, 0};
            }
            break;
        case 2:
            c = { (v - vmin) / dv, 0, (vmax - v) / dv};
            break;
        case 3:
            r = (v - vmin) / dv;
            c = {r, r, r};
            break;
        case 4:
            if (v < (vmin + dv / 6))          {   c = {1, 6 * (v - vmin) / dv, 0};                   }
            else if (v < (vmin + 2 * dv / 6)) {   c = { 1 + 6 * (vmin + dv / 6 - v) / dv, 1, 0 };    }
            else if (v < (vmin + 3 * dv / 6)) {   c = { 0, 1, 6 * (v - vmin - 2 * dv / 6) / dv};     }
            else if (v < (vmin + 4 * dv / 6)) {   c = { 0, 1 + 6 * (vmin + 3 * dv / 6 - v) / dv, 1}; }
            else if (v < (vmin + 5 * dv / 6)) {   c = { 6 * (v - vmin - 4 * dv / 6) / dv, 0, 1};     }
            else                              {   c = { 1, 0, 1 + 6 * (vmin + 5 * dv / 6 - v) / dv };}
            break;
        case 5:
            c = { (v - vmin) / (vmax - vmin), 1, 0 };
            break;
        case 6:
            c.r = (v - vmin) / (vmax - vmin);
            c.g = (vmax - v) / (vmax - vmin);
            c.b = c.r;
            break;
        case 7:
            if (v < (vmin + 0.25 * dv)) {
                c.r = 0;
                c.g = 4 * (v - vmin) / dv;
                c.b = 1 - c.g;
            }
            else if (v < (vmin + 0.5 * dv)) {
                c.r = 4 * (v - vmin - 0.25 * dv) / dv;
                c.g = 1 - c.r;
                c.b = 0;
            }
            else if (v < (vmin + 0.75 * dv)) {
                c.g = 4 * (v - vmin - 0.5 * dv) / dv;
                c.r = 1 - c.g;
                c.b = 0;
            }
            else {
                c.r = 0;
                c.b = 4 * (v - vmin - 0.75 * dv) / dv;
                c.g = 1 - c.b;
            }
            break;
        case 8:
            if (v < (vmin + 0.5 * dv)) {
                c.r = 2 * (v - vmin) / dv;
                c.g = c.r;
                c.b = c.r;
            }
            else {
                c.r = 1 - 2 * (v - vmin - 0.5 * dv) / dv;
                c.g = c.r;
                c.b = c.r;
            }
            break;
        case 9:
            if (v < (vmin + dv / 3)) {
                c.b = 3 * (v - vmin) / dv;
                c.g = 0;
                c.r = 1 - c.b;
            }
            else if (v < (vmin + 2 * dv / 3)) {
                c.r = 0;
                c.g = 3 * (v - vmin - dv / 3) / dv;
                c.b = 1;
            }
            else {
                c.r = 3 * (v - vmin - 2 * dv / 3) / dv;
                c.g = 1 - c.r;
                c.b = 1;
            }
            break;
        case 10:
            if (v < (vmin + 0.2 * dv)) {
                c.r = 0;
                c.g = 5 * (v - vmin) / dv;
                c.b = 1;
            }
            else if (v < (vmin + 0.4 * dv)) {
                c.r = 0;
                c.g = 1;
                c.b = 1 + 5 * (vmin + 0.2 * dv - v) / dv;
            }
            else if (v < (vmin + 0.6 * dv)) {
                c.r = 5 * (v - vmin - 0.4 * dv) / dv;
                c.g = 1;
                c.b = 0;
            }
            else if (v < (vmin + 0.8 * dv)) {
                c.r = 1;
                c.g = 1 - 5 * (v - vmin - 0.6 * dv) / dv;
                c.b = 0;
            }
            else {
                c.r = 1;
                c.g = 5 * (v - vmin - 0.8 * dv) / dv;
                c.b = 5 * (v - vmin - 0.8 * dv) / dv;
            }
            break;
        case 11:
            c1.r = 200 / 255;
            c1.g = 60 / 255;
            c1.b = 0 / 255;
            c2.r = 250 / 255;
            c2.g = 160 / 255;
            c2.b = 110 / 255;
            c.r = (c2.r - c1.r) * (v - vmin) / dv + c1.r;
            c.g = (c2.g - c1.g) * (v - vmin) / dv + c1.g;
            c.b = (c2.b - c1.b) * (v - vmin) / dv + c1.b;
            break;
        case 12:
            c1.r = 55 / 255;
            c1.g = 55 / 255;
            c1.b = 45 / 255;
            /* c2.r = 200 / 255; c2.g =  60 / 255; c2.b =   0 / 255; */
            c2.r = 235 / 255;
            c2.g = 90 / 255;
            c2.b = 30 / 255;
            c3.r = 250 / 255;
            c3.g = 160 / 255;
            c3.b = 110 / 255;
            ratio = 0.4;
            vmid = vmin + ratio * dv;
            if (v < vmid) {
                c.r = (c2.r - c1.r) * (v - vmin) / (ratio * dv) + c1.r;
                c.g = (c2.g - c1.g) * (v - vmin) / (ratio * dv) + c1.g;
                c.b = (c2.b - c1.b) * (v - vmin) / (ratio * dv) + c1.b;
            }
            else {
                c.r = (c3.r - c2.r) * (v - vmid) / ((1 - ratio) * dv) + c2.r;
                c.g = (c3.g - c2.g) * (v - vmid) / ((1 - ratio) * dv) + c2.g;
                c.b = (c3.b - c2.b) * (v - vmid) / ((1 - ratio) * dv) + c2.b;
            }
            break;
        case 13:
            c1.r = 0 / 255;
            c1.g = 255 / 255;
            c1.b = 0 / 255;
            c2.r = 255 / 255;
            c2.g = 150 / 255;
            c2.b = 0 / 255;
            c3.r = 255 / 255;
            c3.g = 250 / 255;
            c3.b = 240 / 255;
            ratio = 0.3;
            vmid = vmin + ratio * dv;
            if (v < vmid) {
                c.r = (c2.r - c1.r) * (v - vmin) / (ratio * dv) + c1.r;
                c.g = (c2.g - c1.g) * (v - vmin) / (ratio * dv) + c1.g;
                c.b = (c2.b - c1.b) * (v - vmin) / (ratio * dv) + c1.b;
            }
            else {
                c.r = (c3.r - c2.r) * (v - vmid) / ((1 - ratio) * dv) + c2.r;
                c.g = (c3.g - c2.g) * (v - vmid) / ((1 - ratio) * dv) + c2.g;
                c.b = (c3.b - c2.b) * (v - vmid) / ((1 - ratio) * dv) + c2.b;
            }
            break;
        case 14:
            c.r = 1;
            c.g = 1 - (v - vmin) / dv;
            c.b = 0;
            break;
        case 15:
            if (v < (vmin + 0.25 * dv)) {
                c.r = 0;
                c.g = 4 * (v - vmin) / dv;
                c.b = 1;
            }
            else if (v < (vmin + 0.5 * dv)) {
                c.r = 0;
                c.g = 1;
                c.b = 1 - 4 * (v - vmin - 0.25 * dv) / dv;
            }
            else if (v < (vmin + 0.75 * dv)) {
                c.r = 4 * (v - vmin - 0.5 * dv) / dv;
                c.g = 1;
                c.b = 0;
            }
            else {
                c.r = 1;
                c.g = 1;
                c.b = 4 * (v - vmin - 0.75 * dv) / dv;
            }
            break;
        case 16:
            if (v < (vmin + 0.5 * dv)) {
                c.r = 0;
                c.g = 2 * (v - vmin) / dv;
                c.b = 1 - 2 * (v - vmin) / dv;
            }
            else {
                c.r = 2 * (v - vmin - 0.5 * dv) / dv;
                c.g = 1 - 2 * (v - vmin - 0.5 * dv) / dv;
                c.b = 0;
            }
            break;
        case 17:
            if (v < (vmin + 0.5 * dv)) {
                c.r = 1;
                c.g = 1 - 2 * (v - vmin) / dv;
                c.b = 2 * (v - vmin) / dv;
            }
            else {
                c.r = 1 - 2 * (v - vmin - 0.5 * dv) / dv;
                c.g = 2 * (v - vmin - 0.5 * dv) / dv;
                c.b = 1;
            }
            break;
        case 18:
            c.r = 0;
            c.g = (v - vmin) / (vmax - vmin);
            c.b = 1;
            break;
        case 19:
            c.r = (v - vmin) / (vmax - vmin);
            c.g = c.r;
            c.b = 1;
            break;
        case 20:
            c1.r = 0 / 255;
            c1.g = 160 / 255;
            c1.b = 0 / 255;
            c2.r = 180 / 255;
            c2.g = 220 / 255;
            c2.b = 0 / 255;
            c3.r = 250 / 255;
            c3.g = 220 / 255;
            c3.b = 170 / 255;
            ratio = 0.3;
            vmid = vmin + ratio * dv;
            if (v < vmid) {
                c.r = (c2.r - c1.r) * (v - vmin) / (ratio * dv) + c1.r;
                c.g = (c2.g - c1.g) * (v - vmin) / (ratio * dv) + c1.g;
                c.b = (c2.b - c1.b) * (v - vmin) / (ratio * dv) + c1.b;
            }
            else {
                c.r = (c3.r - c2.r) * (v - vmid) / ((1 - ratio) * dv) + c2.r;
                c.g = (c3.g - c2.g) * (v - vmid) / ((1 - ratio) * dv) + c2.g;
                c.b = (c3.b - c2.b) * (v - vmid) / ((1 - ratio) * dv) + c2.b;
            }
            break;
        case 21:
            c1.r = 255 / 255;
            c1.g = 255 / 255;
            c1.b = 200 / 255;
            c2.r = 150 / 255;
            c2.g = 150 / 255;
            c2.b = 255 / 255;
            c.r = (c2.r - c1.r) * (v - vmin) / dv + c1.r;
            c.g = (c2.g - c1.g) * (v - vmin) / dv + c1.g;
            c.b = (c2.b - c1.b) * (v - vmin) / dv + c1.b;
            break;
        case 22:
            c.r = 1 - (v - vmin) / dv;
            c.g = 1 - (v - vmin) / dv;
            c.b = (v - vmin) / dv;
            break;
        case 23:
            if (v < (vmin + 0.5 * dv)) {
                c.r = 1;
                c.g = 2 * (v - vmin) / dv;
                c.b = c.g;
            }
            else {
                c.r = 1 - 2 * (v - vmin - 0.5 * dv) / dv;
                c.g = c.r;
                c.b = 1;
            }
            break;
        case 24:
            if (v < (vmin + 0.5 * dv)) {
                c.r = 2 * (v - vmin) / dv;
                c.g = c.r;
                c.b = 1 - c.r;
            }
            else {
                c.r = 1;
                c.g = 1 - 2 * (v - vmin - 0.5 * dv) / dv;
                c.b = 0;
            }
            break;
        case 25:
            if (v < (vmin + dv / 3)) {
                c.r = 0;
                c.g = 3 * (v - vmin) / dv;
                c.b = 1;
            }
            else if (v < (vmin + 2 * dv / 3)) {
                c.r = 3 * (v - vmin - dv / 3) / dv;
                c.g = 1 - c.r;
                c.b = 1;
            }
            else {
                c.r = 1;
                c.g = 0;
                c.b = 1 - 3 * (v - vmin - 2 * dv / 3) / dv;
            }
            break;
    }
    return c;
}