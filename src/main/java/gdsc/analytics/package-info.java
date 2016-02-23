/**
 * <ul>
 * <li>Copyright (c) 2010 Daniel Murphy
 * <li>Copyright (c) 2016 Alex Herbert
 * </ul>
 * 
 * The core of this packgage is based upon JGoogleAnalyticsTracker by Daniel Murphy.
 * A similar package is JGoogleAnalytics by Siddique Hameed.
 * <p>
 * The JGoogleAnalyticsTracker code has been modified to allow better
 * use of the session tracking within Google Analytics using session
 * timeout. The slf4j dependency was removed. 
 * <p>
 * Since the code will only be used within a Java application
 * the referral, search referral and campaign functionality has been
 * removed to simplify the analytics and allow caching most of the
 * constructed analytics URL. The code is redistributed under the
 * original MIT licence.
 * <p>
 * MIT Licence:
 * <p>
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject to
 * the following conditions:
 * <p>
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 * <p>
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 * <p>
 * @see https://code.google.com/archive/p/jgoogleanalyticstracker/
 * @see https://github.com/siddii/jgoogleanalytics
 */
package gdsc.analytics;
