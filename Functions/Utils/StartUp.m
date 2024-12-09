
%% Clean up

sca; clear all; close all; clc; commandwindow;

AssertOpenGL;

%% Set random seed

[seed, whichGen] = ClockRandSeed;

%% Initialize garbage collector????

%% Screen preferences

% Suppress non-critical tests and warnings      (TEST ELSEWHERE!)
Screen('Preference', 'VisualDebuglevel',    3);
Screen('Preference', 'SuppressAllWarnings', 1);
Screen('Preference',  'SkipSyncTests',      2);

% PTB preparation
KbName('UnifyKeyNames');
RestrictKeysForKbCheck([]);
ListenChar(0);
